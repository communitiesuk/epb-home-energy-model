use crate::core::controls::time_control::SetpointOrCombinationControl;
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::heating_systems::instant_elec_heater::InstantElecHeater;
use crate::corpus::TempInternalAirFn;
use crate::hem_core::simulation_time::SimulationTimeIteration;
use anyhow::bail;
use delegate::delegate;

// BS EN 1264 limits temperature to 15C above room temperature in peripheral areas. This is 36C assuming a room temperature of 21C
const MAX_TEMPERATURE: f64 = 36.;

struct DryElectricUnderfloorHeating {
    electric_heater: InstantElecHeater,
    total_emitter_floor_area: f64,
}

impl DryElectricUnderfloorHeating {
    pub fn new(
        rated_power: f64,
        emitter_floor_area: f64,
        frac_convective: f64,
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        zone_area: f64,
        simulation_timestep: f64,
        control: SetpointOrCombinationControl,
        c: Option<f64>,
        c_per_m2: Option<f64>,
        n: f64,
        thermal_mass: Option<f64>,
        thermal_mass_per_m2: Option<f64>,
        initial_temperature: Option<f64>,
    ) -> anyhow::Result<Self> {
        let initial_temperature = initial_temperature.unwrap_or(20.);

        let c = if let Some(c_per_m2) = c_per_m2 {
            Some(c_per_m2 * emitter_floor_area)
        } else {
            c
        };

        let thermal_mass = if let Some(thermal_mass_per_m2) = thermal_mass_per_m2 {
            Some(thermal_mass_per_m2 * emitter_floor_area)
        } else {
            thermal_mass
        };

        if thermal_mass.is_none() {
            bail!("Thermal mass is missing");
        }

        if c.is_none() {
            bail!("Constant from characteristic equation of emitters is missing");
        }

        // check floor area validity
        if emitter_floor_area > zone_area {
            bail!("Total UFH area ({emitter_floor_area}) is bigger than Zone area ({zone_area})");
        }

        Ok(Self {
            electric_heater: InstantElecHeater::new(
                rated_power,
                frac_convective,
                energy_supply_connection,
                temp_internal_air_fn,
                simulation_timestep,
                control.into(),
                c,
                None,
                n.into(),
                thermal_mass,
                None,
                initial_temperature,
                MAX_TEMPERATURE.into(),
            )?,
            total_emitter_floor_area: emitter_floor_area,
        })
    }

    pub fn total_emitter_floor_area(&self) -> f64 {
        self.total_emitter_floor_area
    }

    delegate! {
        to self.electric_heater {
            pub fn demand_energy(&self, energy_demand: f64, simtime: SimulationTimeIteration) -> anyhow::Result<f64>;
            pub fn temp_setpnt(&self, simtime: &SimulationTimeIteration) -> Option<f64>;
            pub fn in_required_period(&self, simtime: &SimulationTimeIteration) -> Option<bool>;
            pub fn frac_convective(&self) -> f64;
            pub fn energy_output_min(&self) -> anyhow::Result<f64>;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
    use crate::hem_core::simulation_time::SimulationTime;
    use crate::input::FuelType;
    use approx::assert_relative_eq;
    use parking_lot::RwLock;
    use rstest::*;
    use std::sync::Arc;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 4., 1.)
    }

    #[fixture]
    fn energy_supply_connection(simulation_time: SimulationTime) -> EnergySupplyConnection {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        EnergySupply::connection(energy_supply, "shower").unwrap()
    }

    #[fixture]
    fn control(simulation_time: SimulationTime) -> SetpointOrCombinationControl {
        SetpointOrCombinationControl::SetpointTime(
            SetpointTimeControl::new(
                vec![Some(21.0), Some(21.0), None, Some(21.0)],
                0,
                1.,
                None,
                None,
                simulation_time.step,
            )
            .into(),
        )
    }

    #[fixture]
    fn zone_area() -> f64 {
        80.
    }

    fn create_temp_internal_air_fn(canned_value: f64) -> TempInternalAirFn {
        Arc::new(move || canned_value)
    }

    #[fixture]
    fn temp_internal_air_fn() -> TempInternalAirFn {
        create_temp_internal_air_fn(20.)
    }

    #[fixture]
    fn heater(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        zone_area: f64,
        control: SetpointOrCombinationControl,
        simulation_time: SimulationTime,
    ) -> DryElectricUnderfloorHeating {
        DryElectricUnderfloorHeating::new(
            50.,
            80.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            zone_area,
            simulation_time.step,
            control,
            Some(1.3),
            None,
            1.2,
            Some(0.14),
            None,
            None,
        )
        .unwrap()
    }

    #[rstest]
    fn test_missing_parameters(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        zone_area: f64,
        control: SetpointOrCombinationControl,
        simulation_time: SimulationTime,
    ) {
        assert!(DryElectricUnderfloorHeating::new(
            50.,
            80.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            zone_area,
            simulation_time.step,
            control,
            None,
            None,
            1.2,
            Some(0.14),
            None,
            None
        )
        .is_err());
    }

    #[rstest]
    /// Test that the constructor throws if the emitter area is above the floor area
    fn test_emitter_area_above_floor_area(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        zone_area: f64,
        control: SetpointOrCombinationControl,
        simulation_time: SimulationTime,
    ) {
        assert!(DryElectricUnderfloorHeating::new(
            50.,
            100.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            zone_area,
            simulation_time.step,
            control,
            Some(0.13),
            None,
            1.2,
            Some(0.14),
            None,
            None
        )
        .is_err());
    }

    #[rstest]
    /// Test that the energy supplied is the same with setting either c or c_per_m2
    fn test_c_per_m2(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        zone_area: f64,
        control: SetpointOrCombinationControl,
        simulation_time: SimulationTime,
    ) {
        let heater1 = DryElectricUnderfloorHeating::new(
            50.,
            80.,
            0.4,
            energy_supply_connection.clone(),
            temp_internal_air_fn.clone(),
            zone_area,
            simulation_time.step,
            control.clone(),
            None,
            Some(1.3 / 80.),
            1.2,
            Some(0.14),
            None,
            None,
        )
        .unwrap();

        let heater2 = DryElectricUnderfloorHeating::new(
            50.,
            80.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            zone_area,
            simulation_time.step,
            control,
            Some(1.3),
            None,
            1.2,
            Some(0.14),
            None,
            None,
        )
        .unwrap();

        let simtime = simulation_time.iter().next().unwrap();

        assert_eq!(
            heater1.demand_energy(20., simtime).unwrap(),
            heater2.demand_energy(20., simtime).unwrap()
        );
    }

    #[rstest]
    /// Test that the energy supplied is the same with setting either thermal_mass or thermal_mass_per_m2
    fn test_thermal_mass_per_m2(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        zone_area: f64,
        control: SetpointOrCombinationControl,
        simulation_time: SimulationTime,
    ) {
        let heater1 = DryElectricUnderfloorHeating::new(
            50.,
            80.,
            0.4,
            energy_supply_connection.clone(),
            temp_internal_air_fn.clone(),
            zone_area,
            simulation_time.step,
            control.clone(),
            Some(1.3),
            None,
            1.2,
            Some(0.14),
            None,
            None,
        )
        .unwrap();

        let heater2 = DryElectricUnderfloorHeating::new(
            50.,
            80.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            zone_area,
            simulation_time.step,
            control,
            Some(1.3),
            None,
            1.2,
            Some(0.14),
            Some(0.14 / 80.),
            None,
        )
        .unwrap();

        let simtime = simulation_time.iter().next().unwrap();

        assert_eq!(
            heater1.demand_energy(20., simtime).unwrap(),
            heater2.demand_energy(20., simtime).unwrap()
        );
    }

    #[rstest]
    /// Test that DryElectricUnderfloorHeating object returns correct energy supplied
    fn test_demand_energy(heater: DryElectricUnderfloorHeating, simulation_time: SimulationTime) {
        let energy_demands = [40.0, 100.0, 30.0, 20.0];
        let expected_energy_supplied = [
            34.66394952823805,
            36.214903433118764,
            29.410138832041746,
            19.112267824658698,
        ];

        for (t_idx, simtime) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                heater
                    .demand_energy(energy_demands[t_idx], simtime)
                    .unwrap(),
                expected_energy_supplied[t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    /// Test that DryElectricUnderfloorHeating throws with no thermal mass
    fn test_demand_energy_no_thermal_mass(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        zone_area: f64,
        control: SetpointOrCombinationControl,
        simulation_time: SimulationTime,
    ) {
        assert!(DryElectricUnderfloorHeating::new(
            50.,
            80.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            zone_area,
            simulation_time.step,
            control,
            Some(1.3),
            None,
            1.2,
            None,
            None,
            None,
        )
        .is_err());
    }

    #[rstest]
    fn test_temp_setpnt(heater: DryElectricUnderfloorHeating, simulation_time: SimulationTime) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                heater.temp_setpnt(&t_it),
                [21.0.into(), 21.0.into(), None, 21.0.into()][t_idx]
            );
        }
    }

    #[rstest]
    fn test_in_required_period(
        heater: DryElectricUnderfloorHeating,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                heater.in_required_period(&t_it),
                [true, true, false, true][t_idx].into()
            );
        }
    }

    #[rstest]
    fn test_frac_convective(heater: DryElectricUnderfloorHeating) {
        assert_eq!(heater.frac_convective(), 0.4);
    }

    #[rstest]
    fn test_energy_output_min(heater: DryElectricUnderfloorHeating) {
        assert_eq!(heater.energy_output_min().unwrap(), 0.0);
    }

    #[rstest]
    /// Test that total_emitter_floor_area returns the correct floor area
    fn test_total_emitter_floor_area(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        zone_area: f64,
        control: SetpointOrCombinationControl,
        simulation_time: SimulationTime,
    ) {
        let heater = DryElectricUnderfloorHeating::new(
            50.,
            80.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            zone_area,
            simulation_time.step,
            control,
            Some(1.3),
            None,
            1.,
            Some(1.),
            None,
            None,
        )
        .unwrap();

        assert_eq!(heater.total_emitter_floor_area(), 80.);
    }
}
