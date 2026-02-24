use crate::compare_floats::{max_of_2, min_of_2};
use crate::external_conditions::ExternalConditions;
use crate::input::BatteryLocation;
use crate::input::ElectricBattery as ElectricBatteryInput;
use crate::simulation_time::SimulationTimeIteration;
use atomic_float::AtomicF64;
use std::sync::atomic::Ordering;
use std::sync::Arc;

/// An object to represent an electric battery system
#[derive(Debug)]
pub struct ElectricBattery {
    /// the maximum capacity of the battery
    _capacity: f64,
    /// charge/discharge round trip efficiency of battery system (between 0 & 1)
    charge_discharge_efficiency: f64,
    minimum_charge_rate: f64,
    maximum_charge_rate: f64,
    maximum_discharge_rate: f64,
    battery_location: BatteryLocation,
    grid_charging_possible: bool,
    simulation_timestep: f64,
    external_conditions: Arc<ExternalConditions>,
    /// the current energy stored in the battery at the
    /// end of hour (kWh)
    current_energy_stored: AtomicF64,
    total_time_charging_current_timestep: AtomicF64,
    state_of_health: f64,
    max_capacity: f64,
    one_way_battery_efficiency: f64,
    reverse_one_way_battery_efficiency: f64,
}

/// Arguments:
/// * `capacity` - the maximum capacity of the battery (kWh)
/// * `charge_discharge_efficiency` - charge/discharge round trip efficiency of battery
///                                       system (between 0 & 1)
/// * `battery_age` - the starting age of the battery in (years)
/// * `minimum_charge_rate` - the minimum charge rate one way trip the battery allows (kW)
/// * `maximum_charge_rate` - the maximum charge rate one way trip the battery allows (kW)
/// * `maximum_discharge_rate` - the maximum discharge rate one way trip the battery allows (kW)
/// * `battery_location` - Location of battery (outside or inside)
/// * `grid_charging_possible` - Is charging from the grid possible?
/// * `simulation_timestep` - timestep of the simulation time
/// * `external_conditions` - reference to ExternalConditions object
impl ElectricBattery {
    pub fn new(
        capacity: f64,
        charge_discharge_efficiency: f64,
        battery_age: f64,
        minimum_charge_rate: f64,
        maximum_charge_rate: f64,
        maximum_discharge_rate: f64,
        battery_location: BatteryLocation,
        grid_charging_possible: bool,
        simulation_timestep: f64,
        external_conditions: Arc<ExternalConditions>,
    ) -> Self {
        let state_of_health = max_of_2(Self::state_of_health_equ(battery_age), 0.);
        // below code, in combination with charge_discharge_efficiency minimum value being 0, is
        // erroneous and leads to NaN error further on (upstream python throws a ZeroDivisionError
        // immediately), reported upstream
        let one_way_battery_efficiency = charge_discharge_efficiency.powf(0.5);
        let reverse_one_way_battery_efficiency = 1. / one_way_battery_efficiency;

        Self {
            _capacity: capacity,
            charge_discharge_efficiency,
            minimum_charge_rate,
            maximum_charge_rate,
            maximum_discharge_rate,
            battery_location,
            grid_charging_possible,
            simulation_timestep,
            external_conditions,
            current_energy_stored: Default::default(),
            total_time_charging_current_timestep: Default::default(),
            state_of_health,
            // Calculate max capacity based on battery original capacity * state of health
            max_capacity: capacity * state_of_health,
            one_way_battery_efficiency,
            reverse_one_way_battery_efficiency,
        }
    }

    pub fn from_input(
        input: &ElectricBatteryInput,
        simulation_timestep: f64,
        external_conditions: Arc<ExternalConditions>,
    ) -> Self {
        let ElectricBatteryInput {
            capacity,
            charge_discharge_efficiency_round_trip,
            battery_age,
            minimum_charge_rate_one_way_trip,
            maximum_charge_rate_one_way_trip,
            maximum_discharge_rate_one_way_trip,
            battery_location,
            grid_charging_possible,
            ..
        } = input;
        Self::new(
            *capacity,
            *charge_discharge_efficiency_round_trip,
            *battery_age,
            *minimum_charge_rate_one_way_trip,
            *maximum_charge_rate_one_way_trip,
            *maximum_discharge_rate_one_way_trip,
            *battery_location,
            *grid_charging_possible,
            simulation_timestep,
            external_conditions,
        )
    }

    /// Equation for charge rate as function of state of charge (used for calculating the maximum charge rate)
    /// We included this as a function in anticipation to the dependency between the charge/discharge rate with the
    /// state of charge. Having considered the evidence available, we don't think we have solid enough basis to
    /// propose an equation to model this dependency yet.
    /// TODO (from Python) Revisit available research to improve charge/discharge rate dependency with state of charge
    ///      Consider paper "Joint State of Charge (SOC) and State of Health (SOH) Estimation for Lithium-Ion Batteries
    ///      Packs of Electric Vehicles. Principal Author: Panpan Hu
    fn charge_rate_soc_equ(_state_of_charge: f64) -> f64 {
        1.
    }

    /// Equation for discharge rate as function of state of charge (used for calculating the maximum discharge rate)
    /// We included this as a function in anticipation to the dependency between the charge/discharge rate with the
    /// state of charge. Having considered the evidence available, we don't think we have solid enough basis to
    /// propose an equation to model this dependency yet.
    /// TODO (from Python) Revisit available research to improve charge/discharge rate dependency with state of charge
    ///      Consider paper "Joint State of Charge (SOC) and State of Health (SOH) Estimation for Lithium-Ion Batteries
    ///      Packs of Electric Vehicles. Principal Author: Panpan Hu
    fn discharge_rate_soc_equ(_state_of_charge: f64) -> f64 {
        1.
    }

    /// Equation for battery capacity as function of external air temperature (only used if battery is outside)
    /// Based on manufacturer data (graph): https://www.bonnenbatteries.com/the-effect-of-low-temperature-on-lithium-batteries/
    /// TODO (from Python) Revisit available research to improve capacity temperature dependency
    fn capacity_temp_equ(air_temp: f64) -> f64 {
        if air_temp > 20. {
            1.
        } else {
            0.8496 + 0.01208 * air_temp - 0.000228 * air_temp.powi(2)
        }
    }

    /// Equation for battery state of health as function of battery age
    /// Based on manufacturer guarantee of 60% remaining capacity after 10 years - i.e. 4% drop per year.
    /// Guaranteed performance is probably quite conservative. Seek to refine this in future.
    /// TODO (from Python) Seek less conservative source
    fn state_of_health_equ(battery_age: f64) -> f64 {
        -0.04 * battery_age + 1.
    }

    pub(crate) fn get_charge_efficiency(&self, simtime: SimulationTimeIteration) -> f64 {
        let air_temp_capacity_factor = self.limit_capacity_due_to_temp(simtime);

        self.one_way_battery_efficiency * self.state_of_health * air_temp_capacity_factor
    }

    pub(crate) fn get_discharge_efficiency(&self, simtime: SimulationTimeIteration) -> f64 {
        let air_temp_capacity_factor = self.limit_capacity_due_to_temp(simtime);

        self.reverse_one_way_battery_efficiency * self.state_of_health * air_temp_capacity_factor
    }

    pub(crate) fn is_grid_charging_possible(&self) -> bool {
        self.grid_charging_possible
    }

    pub(crate) fn get_state_of_charge(&self) -> f64 {
        self.current_energy_stored.load(Ordering::SeqCst) / self.max_capacity
    }

    pub(crate) fn get_max_capacity(&self) -> f64 {
        self.max_capacity
    }

    /// Arguments:
    /// energy_flow -- the supply (-ve) or demand (+ve) to/on the electric battery (kWh)
    /// charging_from_grid        -- Charging from charging_from_grid, not PV.
    /// Other variables:
    /// energy_available_to_charge_battery -- the total energy that could charge/discharge the battery
    /// including losses from charging efficiency (kWh)
    /// current_energy_stored_unconstrained -- current energy stored in battery + the total energy that
    /// would charge/discharge the battery without minimum/maximum
    /// constraints of the battery (kWh)
    /// energy_accepted_by_battery -- the total energy the battery is able to supply or charge (kWh)
    pub(crate) fn charge_discharge_battery(
        &self,
        energy_flow: f64,
        charging_from_grid: bool,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        let timestep = self.simulation_timestep;
        let total_time_charging_current_timestep = self
            .total_time_charging_current_timestep
            .load(Ordering::SeqCst);

        if (timestep < total_time_charging_current_timestep
            || is_close!(
                timestep,
                total_time_charging_current_timestep,
                rel_tol = 1e-09
            ))
            && energy_flow < 0.
        {
            // No more scope for charging
            return 0.;
        }

        let mut max_charge_rate: Option<f64> = None;

        // Calculate State of Charge (SoC)
        let state_of_charge = Self::get_state_of_charge(self);

        // Calculate the impact on the battery capacity of air temperature
        let air_temp_capacity_factor = self.limit_capacity_due_to_temp(simtime);

        // Ensure state_of_charge is between 0 and 1
        let state_of_charge = min_of_2(state_of_charge, 1.);
        let state_of_charge = max_of_2(state_of_charge, 0.);

        let energy_available_to_charge_battery = if energy_flow < 0. {
            // Charging battery
            // Convert energy_flow (in kWh) to a power (in kW) by dividing energy by time available for charging (in hours)
            let energy_flow_power = energy_flow
                / (timestep
                    - self
                        .total_time_charging_current_timestep
                        .load(Ordering::SeqCst));
            // If supply is less than minimum charge rate, do not add charge to the battery
            if !charging_from_grid && -energy_flow_power < self.minimum_charge_rate {
                0.
            } else {
                let (max_charge, max_charge_rate_temp) = self.calculate_max_charge(state_of_charge);
                max_charge_rate = max_charge_rate_temp;
                min_of_2(
                    -energy_flow
                        * self.one_way_battery_efficiency
                        * self.state_of_health
                        * air_temp_capacity_factor,
                    max_charge,
                )
            }
        } else {
            // Discharging battery
            let max_discharge = if charging_from_grid {
                0.
            } else {
                // max charge energy the battery can supply in the timestep
                self.calculate_max_discharge(state_of_charge)
            };

            max_of_2(-energy_flow, max_discharge) / self.one_way_battery_efficiency
                * self.state_of_health
                * air_temp_capacity_factor // reductions due to state of health (i.e age) and cold temperature applied here
        };

        // Charge/discharge the battery by the amount available
        let current_energy_stored_unconstrained =
            self.current_energy_stored.load(Ordering::SeqCst) + energy_available_to_charge_battery;
        let prev_energy_stored = self.current_energy_stored.load(Ordering::SeqCst);
        // Energy stored cannot be > battery capacity or < 0
        self.current_energy_stored.store(
            min_of_2(
                self.max_capacity,
                max_of_2(0., current_energy_stored_unconstrained),
            ),
            Ordering::SeqCst,
        );
        let energy_accepted_by_battery =
            self.current_energy_stored.load(Ordering::SeqCst) - prev_energy_stored;

        // Return the supply/demand energy the battery can accept (including charging/discharging losses)
        if energy_flow < 0. {
            // Charging battery
            // Calculate charging time of battery
            let time_charging_current_load = match max_charge_rate {
                Some(max_charge_rate) if max_charge_rate >= 0. => min_of_2(
                    energy_accepted_by_battery
                        / max_charge_rate
                        / self.one_way_battery_efficiency
                        / self.state_of_health
                        / air_temp_capacity_factor,
                    timestep
                        - self
                            .total_time_charging_current_timestep
                            .load(Ordering::SeqCst),
                ),
                // TODO (from Python) some of these adjustment factors are probably adjusting the energy_accepted_by_battery while others are adjusting max_charge_rate,
                //  and I think it would clearer to calculate adjusted versions of those two variables separately.
                _ => 0.,
            };

            self.total_time_charging_current_timestep
                .fetch_add(time_charging_current_load, Ordering::SeqCst);

            -energy_accepted_by_battery / self.one_way_battery_efficiency
        } else {
            // Discharging battery
            -energy_accepted_by_battery * self.one_way_battery_efficiency
        }
    }

    fn limit_capacity_due_to_temp(&self, simtime: SimulationTimeIteration) -> f64 {
        let air_temp = match self.battery_location {
            BatteryLocation::Outside => self.external_conditions.air_temp(&simtime),
            BatteryLocation::Inside => 20., // Fixed for now, but if we add zone location of battery (if inside) as an input, we could look this up in future.
        };

        max_of_2(0., min_of_2(1., Self::capacity_temp_equ(air_temp)))
    }

    /// Calculate the maximum rate of charge rate based on the current state of charge of the battery
    /// Arguments:
    /// * `state_of_charge` - the current state of charge (0 - 1) of the battery
    fn calculate_max_charge(&self, state_of_charge: f64) -> (f64, Option<f64>) {
        let charge_factor_for_soc = Self::charge_rate_soc_equ(state_of_charge);
        let max_charge_rate = self.maximum_charge_rate * charge_factor_for_soc;
        let max_charge =
            self.maximum_charge_rate * charge_factor_for_soc * self.simulation_timestep;

        (max_charge, Some(max_charge_rate))
    }

    /// Calculate the maximum amount of discharge possible based on the current state of charge of the battery
    /// Arguments:
    /// * `state_of_charge` - the current state of charge (0 - 1) of the battery
    pub(crate) fn calculate_max_discharge(&self, state_of_charge: f64) -> f64 {
        let discharge_factor_for_soc = Self::discharge_rate_soc_equ(state_of_charge);

        -(self.maximum_discharge_rate * discharge_factor_for_soc * self.simulation_timestep)
    }

    /// Calculations to be done at the end of each timestep
    pub(crate) fn timestep_end(&self) {
        // variables below need to be reset at the end of each timestep
        self.total_time_charging_current_timestep
            .store(Default::default(), Ordering::SeqCst);
    }

    pub(crate) fn get_charge_discharge_efficiency(&self) -> f64 {
        self.charge_discharge_efficiency
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::units::Orientation360;
    use crate::external_conditions::{
        DaylightSavingsConfig, ShadingObject, ShadingObjectType, ShadingSegment,
    };
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use rstest::*;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0],
            vec![3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1],
            vec![0., 20., 40., 60., 0., 20., 40., 60.]
                .into_iter()
                .map(Into::into)
                .collect(),
            vec![11., 25., 42., 52., 60., 44., 28., 15.],
            vec![11., 25., 42., 52., 60., 44., 28., 15.],
            vec![0.2; 8],
            51.42,
            -0.75,
            0,
            0,
            Some(0),
            1.,
            Some(1),
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            vec![
                ShadingSegment {
                    start360: Orientation360::create_from_180(180.).unwrap(),
                    end360: Orientation360::create_from_180(135.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(135.).unwrap(),
                    end360: Orientation360::create_from_180(90.).unwrap(),
                    shading_objects: vec![ShadingObject {
                        object_type: ShadingObjectType::Overhang,
                        height: 2.2,
                        distance: 6.,
                    }],
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(90.).unwrap(),
                    end360: Orientation360::create_from_180(45.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(45.).unwrap(),
                    end360: Orientation360::create_from_180(0.).unwrap(),
                    shading_objects: vec![
                        ShadingObject {
                            object_type: ShadingObjectType::Obstacle,
                            height: 40.,
                            distance: 4.,
                        },
                        ShadingObject {
                            object_type: ShadingObjectType::Overhang,
                            height: 3.,
                            distance: 7.,
                        },
                    ],
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(0.).unwrap(),
                    end360: Orientation360::create_from_180(-45.).unwrap(),
                    shading_objects: vec![ShadingObject {
                        object_type: ShadingObjectType::Obstacle,
                        height: 3.,
                        distance: 8.,
                    }],
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(-45.).unwrap(),
                    end360: Orientation360::create_from_180(-90.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(-90.).unwrap(),
                    end360: Orientation360::create_from_180(-135.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(-135.).unwrap(),
                    end360: Orientation360::create_from_180(-180.).unwrap(),
                    ..Default::default()
                },
            ]
            .into(),
        )
    }

    #[fixture]
    fn electric_battery(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) -> ElectricBattery {
        ElectricBattery::new(
            2.,
            0.8,
            3.,
            0.001,
            1.5,
            1.5,
            BatteryLocation::Outside,
            false,
            simulation_time.step,
            Arc::new(external_conditions),
        )
    }

    #[rstest]
    fn test_get_charge_efficiency(
        electric_battery: ElectricBattery,
        simulation_time: SimulationTime,
    ) {
        let simulation_time = simulation_time.iter().next().unwrap();

        assert_relative_eq!(
            electric_battery.get_charge_efficiency(simulation_time),
            0.6687167004967052
        );
    }

    #[rstest]
    fn test_get_discharge_efficiency(
        electric_battery: ElectricBattery,
        simulation_time: SimulationTime,
    ) {
        let simulation_time = simulation_time.iter().next().unwrap();

        assert_relative_eq!(
            electric_battery.get_discharge_efficiency(simulation_time),
            0.8358958756208815
        );
    }

    #[rstest]
    fn test_get_charge_discharge_efficiency(electric_battery: ElectricBattery) {
        assert_relative_eq!(electric_battery.get_charge_discharge_efficiency(), 0.8);
    }

    #[rstest]
    fn test_get_max_capacity(electric_battery: ElectricBattery) {
        assert_relative_eq!(electric_battery.get_max_capacity(), 1.76);
    }

    #[rstest]
    fn test_timestep_end(electric_battery: ElectricBattery) {
        electric_battery
            .total_time_charging_current_timestep
            .store(10., Ordering::SeqCst);
        electric_battery.timestep_end();

        assert_eq!(
            electric_battery
                .total_time_charging_current_timestep
                .load(Ordering::SeqCst),
            0.
        );
    }

    #[rstest]
    fn test_limit_capacity_due_to_temp(
        electric_battery: ElectricBattery,
        simulation_time: SimulationTime,
    ) {
        let simulation_time = simulation_time.iter().next().unwrap();

        assert_relative_eq!(
            electric_battery.limit_capacity_due_to_temp(simulation_time),
            0.8496
        );

        let mut electric_battery_inside = electric_battery;
        electric_battery_inside.battery_location = BatteryLocation::Inside;

        assert_relative_eq!(
            electric_battery_inside.limit_capacity_due_to_temp(simulation_time),
            1.
        );
    }

    #[rstest]
    fn test_capacity_temp_equ() {
        assert_relative_eq!(ElectricBattery::capacity_temp_equ(5.), 0.9043);
        assert_relative_eq!(ElectricBattery::capacity_temp_equ(10.), 0.9476);
        assert_relative_eq!(ElectricBattery::capacity_temp_equ(20.), 1.);
        assert_relative_eq!(ElectricBattery::capacity_temp_equ(30.), 1.);
    }

    #[rstest]
    fn test_charge_discharge_battery(
        electric_battery: ElectricBattery,
        simulation_time: SimulationTime,
    ) {
        let simulation_time = simulation_time.iter().next().unwrap();
        // supply to battery exceeds limit
        assert_relative_eq!(
            electric_battery.charge_discharge_battery(-1_000., false, simulation_time),
            -1.6770509831248424,
            max_relative = 1e-7
        );
        // demand on battery exceeds limit
        assert_relative_eq!(
            electric_battery.charge_discharge_battery(1_000., false, simulation_time),
            1.121472,
            max_relative = 1e-7
        );
        // normal charge
        assert_relative_eq!(
            electric_battery.charge_discharge_battery(-0.2, false, simulation_time),
            0.,
            max_relative = 1e-7
        );
        // normal discharge
        assert_relative_eq!(
            electric_battery.charge_discharge_battery(0.1, false, simulation_time),
            0.0747648,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_charge_discharge_battery_below_minimum_charge_rate(
        electric_battery: ElectricBattery,
        simulation_time: SimulationTime,
    ) {
        // Test that the battery is not charged when below the minimum charge rate
        let simulation_time = simulation_time.iter().next().unwrap();
        let mut electric_battery = electric_battery;
        electric_battery.minimum_charge_rate = 20.;
        electric_battery.charge_discharge_battery(-10., false, simulation_time);

        assert_eq!(
            electric_battery.charge_discharge_battery(10., false, simulation_time),
            0.
        );
    }

    #[rstest]
    fn test_charge_discharge_battery_no_grid(
        electric_battery: ElectricBattery,
        simulation_time: SimulationTime,
    ) {
        let simulation_time = simulation_time.iter().next().unwrap();
        electric_battery.charge_discharge_battery(-10., true, simulation_time);

        assert_relative_eq!(
            electric_battery.charge_discharge_battery(10., true, simulation_time),
            0.,
        );
    }

    #[rstest]
    fn test_charge_discharge_battery_grid(
        electric_battery: ElectricBattery,
        simulation_time: SimulationTime,
    ) {
        let simulation_time = simulation_time.iter().next().unwrap();
        electric_battery.charge_discharge_battery(-10., false, simulation_time);

        assert_relative_eq!(
            electric_battery.charge_discharge_battery(10., false, simulation_time),
            1.121472,
        );
    }

    #[rstest]
    fn test_charge_discharge_battery_total_time(
        electric_battery: ElectricBattery,
        simulation_time: SimulationTime,
    ) {
        let simulation_time = simulation_time.iter().next().unwrap();
        electric_battery.charge_discharge_battery(1., false, simulation_time);

        assert_relative_eq!(
            electric_battery
                .total_time_charging_current_timestep
                .load(Ordering::SeqCst),
            0.,
        );

        electric_battery.charge_discharge_battery(-1., false, simulation_time);

        assert_relative_eq!(
            electric_battery
                .total_time_charging_current_timestep
                .load(Ordering::SeqCst),
            2. / 3.
        );

        electric_battery.charge_discharge_battery(1., false, simulation_time);

        assert_relative_eq!(
            electric_battery
                .total_time_charging_current_timestep
                .load(Ordering::SeqCst),
            2. / 3.
        );
    }

    #[rstest]
    fn test_charge_discharge_battery_total_time_no_charge(
        electric_battery: ElectricBattery,
        simulation_time: SimulationTime,
    ) {
        let simulation_time = simulation_time.iter().next().unwrap();
        let mut electric_battery = electric_battery;
        electric_battery.minimum_charge_rate = 20.;
        electric_battery.charge_discharge_battery(1., false, simulation_time);

        assert_relative_eq!(
            electric_battery
                .total_time_charging_current_timestep
                .load(Ordering::SeqCst),
            0.
        );

        electric_battery.charge_discharge_battery(-1., false, simulation_time);

        assert_relative_eq!(
            electric_battery
                .total_time_charging_current_timestep
                .load(Ordering::SeqCst),
            0.
        );

        electric_battery.charge_discharge_battery(1., false, simulation_time);

        assert_relative_eq!(
            electric_battery
                .total_time_charging_current_timestep
                .load(Ordering::SeqCst),
            0.
        );
    }
}
