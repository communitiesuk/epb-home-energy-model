use crate::core::controls::time_control::{per_control, ControlBehaviour};
use crate::core::heating_systems::heat_battery_drycore::{
    convert_to_kwh, HeatBatteryDryCoreCommonBehaviour,
};
use crate::core::heating_systems::heat_battery_drycore::{HeatStorageDryCore, OutputMode};
use crate::{
    core::{controls::time_control::Control, energy_supply::energy_supply::EnergySupplyConnection},
    external_conditions::ExternalConditions,
    input::{ControlLogicType, ElectricStorageHeaterAirFlowType},
    simulation_time::{SimulationTimeIteration, SimulationTimeIterator},
};
use anyhow::bail;
use derivative::Derivative;
use nalgebra::{Vector1, Vector3};
use parking_lot::RwLock;
use std::sync::Arc;

type State = Vector1<f64>;
type Time = f64;

type EnergyOutputState = Vector3<f64>;

// replicates numpys's linspace function
pub(super) fn linspace(start: f64, end: f64, num: i32) -> Vec<f64> {
    let step = (end - start) / f64::from(num - 1);
    (0..num).map(|n| start + (f64::from(n) * step)).collect()
}

#[derive(Derivative)]
#[derivative(Debug)]
pub(crate) struct ElecStorageHeater {
    storage: Arc<RwLock<HeatStorageDryCore>>,
    pwr_instant: f64,
    air_flow_type: ElectricStorageHeaterAirFlowType,
    frac_convective: f64,
    energy_supply_conn: EnergySupplyConnection,
    control: Arc<Control>,
    fan_pwr: f64,
    external_conditions: Arc<ExternalConditions>,
    temp_air: f64,
    zone_setpoint_init: f64,
    #[derivative(Debug = "ignore")]
    zone_internal_air_func: Arc<dyn Fn() -> f64 + Send + Sync>,
    current_energy_profile: RwLock<CurrentEnergyProfile>,
    esh_detailed_results: Option<Arc<RwLock<Vec<StorageHeaterDetailedResult>>>>,
}

#[derive(Clone, Copy, Debug, Default)]
/// A struct to encapsulate energy values for a current step.
struct CurrentEnergyProfile {
    energy_for_fan: f64,
    energy_instant: f64,
    energy_charged: f64,
    energy_delivered: f64,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct StorageHeaterDetailedResult {
    timestep_idx: usize,
    n_units: u32,
    energy_demand: f64,
    energy_delivered: f64,
    energy_instant: f64,
    energy_charged: f64,
    energy_for_fan: f64,
    state_of_charge: f64,
    final_soc: f64,
    time_used_max: f64,
}

impl StorageHeaterDetailedResult {
    pub(crate) fn as_string_values(&self) -> Vec<String> {
        vec![
            self.timestep_idx.to_string(),
            self.n_units.to_string(),
            self.energy_demand.to_string(),
            self.energy_delivered.to_string(),
            self.energy_instant.to_string(),
            self.energy_charged.to_string(),
            self.energy_for_fan.to_string(),
            self.state_of_charge.to_string(),
            self.final_soc.to_string(),
            self.time_used_max.to_string(),
        ]
    }
}

impl ElecStorageHeater {
    /// Arguments:
    /// * `pwr_in` - in kW (Charging)
    /// * `rated_power_instant` - in kW (Instant backup)
    /// * `storage_capacity` - in kWh
    /// * `air_flow_type` - str enum specifying type of Electric Storage Heater:
    ///     * AirFlowType.FAN_ASSISTED
    ///     * AirFlowType.DAMPER_ONLY
    /// * `frac_convective`      - convective fraction for heating (TODO: Check if necessary)
    /// * `fan_pwr`              - Fan power [W]
    /// * `n_units`              - number of units install in zone
    /// * `zone_internal_air_func`  - function that provides access to the internal air temp of the zone
    /// * `energy_supply_conn`   - reference to EnergySupplyConnection object
    /// * `simulation_time`      - reference to SimulationTime object
    /// * `control`              - reference to a control object which must implement is_on() and setpnt() funcs
    /// * `charge_control`       - reference to a ChargeControl object which must implement different logic types
    ///                         for charging the Electric Storage Heaters.
    /// * `dry_core_min_output`       - Data from test showing the output from the storage heater when not actively
    ///                         outputting heat, i.e. case losses only (with units kW)
    /// * `dry_core_max_output`       - Data from test showing the output from the storage heater when it is actively
    ///                         outputting heat, e.g. damper open / fan running (with units kW)
    /// * `external_conditions`  - reference to ExternalConditions object
    pub(crate) fn new(
        pwr_in: f64,
        rated_power_instant: f64,
        storage_capacity: f64,
        air_flow_type: ElectricStorageHeaterAirFlowType,
        frac_convective: f64,
        fan_pwr: f64,
        n_units: u32,
        zone_setpoint_init: f64,
        zone_internal_air_func: Arc<dyn Fn() -> f64 + Send + Sync>,
        energy_supply_conn: EnergySupplyConnection,
        simulation_time: &SimulationTimeIterator,
        control: Arc<Control>,
        charge_control: Arc<Control>,
        dry_core_min_output: Vec<[f64; 2]>,
        dry_core_max_output: Vec<[f64; 2]>,
        external_conditions: Arc<ExternalConditions>,
        state_of_charge_init: f64,
        output_detailed_results: Option<bool>,
    ) -> anyhow::Result<Arc<Self>> {
        let output_detailed_results = output_detailed_results.unwrap_or(false);

        match charge_control.as_ref() {
            Control::Charge(charge) => {
                if charge.logic_type() == ControlLogicType::HeatBattery {
                    bail!("Control logic type HeatBattery is not valid for ElecStorageHeater.")
                }
            }
            _ => bail!("charge_control must be a ChargeControl"),
        }

        let temp_air = zone_internal_air_func();

        let storage = Arc::new(RwLock::new(HeatStorageDryCore::new(
            pwr_in,
            storage_capacity,
            n_units,
            charge_control.clone(),
            dry_core_min_output,
            dry_core_max_output,
            state_of_charge_init,
        )?));

        let heater = Self {
            storage: storage.clone(),
            pwr_instant: rated_power_instant,
            air_flow_type,
            frac_convective,
            energy_supply_conn,
            control,
            fan_pwr,
            external_conditions,
            temp_air,
            zone_setpoint_init,
            zone_internal_air_func,
            current_energy_profile: Default::default(),
            esh_detailed_results: output_detailed_results.then(|| {
                Arc::new(RwLock::new(Vec::with_capacity(
                    simulation_time.total_steps(),
                )))
            }),
        };

        let heater = Arc::new(heater);
        storage.write().set_owner(heater.clone());

        Ok(heater)
    }

    pub(crate) fn temp_setpnt(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<f64> {
        per_control!(self.control.as_ref(), ctrl => { ctrl.setpnt(simulation_time_iteration) })
    }

    pub(crate) fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        per_control!(self.control.as_ref(), ctrl => { ctrl.in_required_period(simulation_time_iteration) })
    }

    pub(crate) fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    /// Calculates the minimum energy that must be delivered based on dry_core_min_output.
    /// :return: np.float64 (minimum energy deliverable in kWh * __n_units).
    pub(crate) fn energy_output_min(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        Ok(self
            .storage
            .read()
            .energy_output(OutputMode::Min, None, None, simulation_time_iteration)?
            .0
            * self.storage.read().n_units() as f64)
    }

    pub(crate) fn energy_output_max(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, f64)> {
        // Calculates the maximum energy that can be delivered based on ESH_max_output.
        // :return: Tuple containing (maximum energy deliverable in kWh, time used in hours).
        self.storage
            .read()
            .energy_output(OutputMode::Max, None, None, simulation_time_iteration)
    }

    pub(crate) fn demand_energy(
        &self,
        energy_demand: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Determines the amount of energy to release based on energy demand, while also handling the
        // energy charging and logging fan energy.
        // :param energy_demand: Energy demand in kWh.
        // :return: Total net energy delivered (including instant heating and fan energy).
        let mut current_profile = self.current_energy_profile.write();

        let timestep = simulation_time_iteration.timestep;
        let n_units: u32 = self.storage.read().n_units();
        let energy_demand = energy_demand / f64::from(n_units);
        current_profile.energy_instant = 0.;

        // Initialize time_used_max and energy_charged_max to default values
        let mut time_used_max = 0.;

        // Calculate minimum energy that can be delivered
        let (q_released_min, _, energy_charged, mut final_soc) = self
            .storage
            .read()
            .energy_output(OutputMode::Min, None, None, simulation_time_iteration)?;
        current_profile.energy_charged = energy_charged;

        let mut q_released_max: Option<f64> = None;

        if q_released_min > energy_demand {
            // Deliver at least the minimum energy
            current_profile.energy_delivered = q_released_min;
            self.storage.write().set_demand_met(q_released_min);
            self.storage.write().set_demand_unmet(0.);
        } else {
            // Calculate maximum energy that can be delivered
            let (q_released_max_value, time_used_max_tmp, energy_charged_max, final_soc_override) =
                self.energy_output_max(simulation_time_iteration)?;
            final_soc = final_soc_override;

            q_released_max = Some(q_released_max_value);
            time_used_max = time_used_max_tmp;

            if q_released_max_value < energy_demand {
                // Deliver as much as possible up to the maximum energy
                current_profile.energy_delivered = q_released_max_value;
                self.storage.write().set_demand_met(q_released_max_value);
                self.storage
                    .write()
                    .set_demand_unmet(energy_demand - q_released_max_value);
                current_profile.energy_charged = energy_charged_max;

                // For now, we assume demand not met from storage is topped-up by
                // the direct top-up heater (if applicable). If still some unmet,
                // this is reported as unmet demand.
                if self.pwr_instant != 0. {
                    current_profile.energy_instant = self
                        .storage
                        .read()
                        .demand_unmet()
                        .min(self.pwr_instant * timestep); // kWh
                    let time_instant = current_profile.energy_instant / self.pwr_instant;
                    time_used_max += time_instant;
                    time_used_max = time_used_max.min(timestep);
                }
            } else {
                // IMPROVED ACCURACY: Use exact energy target instead of linear proration
                // Deliver exactly the demanded energy using the differential equation solver
                let (energy_delivered, time_used_max_tmp, energy_charged, final_soc_override) =
                    self.storage.read().energy_output(
                        OutputMode::Max,
                        None,
                        Some(energy_demand),
                        simulation_time_iteration,
                    )?;

                time_used_max = time_used_max_tmp;
                current_profile.energy_charged = energy_charged;
                final_soc = final_soc_override;

                // The solver should have delivered exactly what we requested
                // (within numerical tolerances)
                current_profile.energy_delivered = energy_delivered.min(energy_demand);
                self.storage
                    .write()
                    .set_demand_met(current_profile.energy_delivered);
                self.storage.write().set_demand_unmet(
                    0.0_f64.max(energy_demand - current_profile.energy_delivered),
                );
            }
        }

        // Ensure energy_delivered does not exceed q_released_max
        let max = q_released_max.unwrap_or(q_released_min);
        current_profile.energy_delivered = current_profile.energy_delivered.min(max);

        // Update state of charge using the final SOC from the differential equation solver
        // The ODE has already integrated the charging and discharging accurately
        self.storage.write().set_state_of_charge(final_soc);

        // Calculate fan energy
        current_profile.energy_for_fan = 0.;

        if self.air_flow_type == ElectricStorageHeaterAirFlowType::FanAssisted
            && q_released_max.is_some()
        {
            let power_for_fan = self.fan_pwr;
            current_profile.energy_for_fan = convert_to_kwh(power_for_fan, time_used_max);
        }

        // Log the energy charged, fan energy, and total energy delivered
        let amount_demanded = f64::from(n_units)
            * (current_profile.energy_charged
                + current_profile.energy_instant
                + current_profile.energy_for_fan);
        self.energy_supply_conn
            .demand_energy(amount_demanded, simulation_time_iteration.index)?;

        // If detailed results flag is set populate with values
        if let Some(esh_detailed_results) = &self.esh_detailed_results {
            let CurrentEnergyProfile {
                energy_for_fan,
                energy_instant,
                energy_charged,
                energy_delivered,
            } = *current_profile;
            let result = StorageHeaterDetailedResult {
                timestep_idx: simulation_time_iteration.index,
                n_units,
                energy_delivered,
                energy_demand,
                energy_instant,
                energy_charged,
                energy_for_fan,
                state_of_charge: self.storage.read().state_of_charge(),
                final_soc,
                time_used_max,
            };
            esh_detailed_results
                .write()
                .insert(simulation_time_iteration.index, result);
        }

        // Return total net energy delivered (discharged + instant heat + fan energy)
        Ok(
            f64::from(n_units)
                * (current_profile.energy_delivered + current_profile.energy_instant),
        )
    }

    pub(crate) fn target_electric_charge(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let storage = self.storage.read();
        storage.target_electric_charge(simulation_time_iteration)
    }

    pub(crate) fn output_esh_results(&self) -> Option<Vec<StorageHeaterDetailedResult>> {
        self.esh_detailed_results
            .as_ref()
            .map(|results| (*results.read()).clone())
    }

    #[cfg(test)]
    fn energy_for_fan(&self) -> f64 {
        self.current_energy_profile.read().energy_for_fan
    }

    #[cfg(test)]
    fn energy_instant(&self) -> f64 {
        self.current_energy_profile.read().energy_instant
    }

    #[cfg(test)]
    fn energy_charged(&self) -> f64 {
        self.current_energy_profile.read().energy_charged
    }

    #[cfg(test)]
    fn energy_delivered(&self) -> f64 {
        self.current_energy_profile.read().energy_delivered
    }
}

impl HeatBatteryDryCoreCommonBehaviour for ElecStorageHeater {
    /// Get temperature for charge control calculations.
    fn get_temp_for_charge_control(&self) -> Option<f64> {
        Some((self.zone_internal_air_func)())
    }

    /// Get zone setpoint for HHRSH calculations.
    fn get_zone_setpoint(&self) -> f64 {
        self.zone_setpoint_init
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::excessive_precision)]
    use super::*;
    use crate::{
        core::{
            controls::time_control::{ChargeControl, SetpointTimeControl},
            energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder},
        },
        external_conditions::{DaylightSavingsConfig, ExternalConditions},
        input::{ControlLogicType, ExternalSensor, FuelType},
        simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator},
    };
    use approx::assert_relative_eq;
    use parking_lot::RwLock;
    use rstest::{fixture, rstest};
    use serde_json::json;

    const EIGHT_DECIMAL_PLACES: f64 = 1e-7;
    const DRY_CORE_MIN_OUTPUT: [[f64; 2]; 3] = [[0.0, 0.0], [0.5, 0.02], [1.0, 0.05]];
    const DRY_CORE_MAX_OUTPUT: [[f64; 2]; 3] = [[0.0, 0.0], [0.5, 1.5], [1.0, 3.0]];

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 24., 1.)
    }

    #[fixture]
    fn simulation_time_iterator(simulation_time: SimulationTime) -> SimulationTimeIterator {
        simulation_time.iter()
    }

    #[fixture]
    fn simulation_time_iteration(
        simulation_time_iterator: SimulationTimeIterator,
    ) -> SimulationTimeIteration {
        simulation_time_iterator.current_iteration()
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> Arc<ExternalConditions> {
        Arc::new(ExternalConditions::new(
            &simulation_time.iter(),
            vec![
                19.0, 0.0, 1.0, 2.0, 5.0, 7.0, 6.0, 12.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
                19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
            ],
            vec![
                3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
                3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
            ],
            vec![
                300., 250., 220., 180., 150., 120., 100., 80., 60., 40., 20., 10., 50., 100., 140.,
                190., 200., 320., 330., 340., 350., 355., 315., 5.,
            ],
            vec![
                0., 0., 0., 0., 35., 73., 139., 244., 320., 361., 369., 348., 318., 249., 225.,
                198., 121., 68., 19., 0., 0., 0., 0., 0.,
            ],
            vec![
                0., 0., 0., 0., 0., 0., 7., 53., 63., 164., 339., 242., 315., 577., 385., 285.,
                332., 126., 7., 0., 0., 0., 0., 0.,
            ],
            vec![
                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
            ],
            51.383,
            -0.783,
            0,
            0,
            Some(0),
            1.,
            Some(1),
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            serde_json::from_value(json!(
                [
                    {"start360": 180, "end360": 135},
                    {"start360": 135, "end360": 90},
                    {"start360": 90, "end360": 45},
                    {"start360": 45, "end360": 0,
                        "shading": [
                            {"type": "obstacle", "height": 10.5, "distance": 12}
                        ]
                    },
                    {"start360": 0, "end360": -45},
                    {"start360": -45, "end360": -90},
                    {"start360": -90, "end360": -135},
                    {"start360": -135, "end360": -180}
                ]
            ))
            .unwrap(),
        ))
    }

    #[fixture]
    fn external_sensor() -> ExternalSensor {
        serde_json::from_value(json!({
            "correlation": [
                {"temperature": 0.0, "max_charge": 1.0},
                {"temperature": 10.0, "max_charge": 0.9},
                {"temperature": 18.0, "max_charge": 0.5}
            ]
        }))
        .unwrap()
    }

    #[fixture]
    fn charge_control_schedule() -> Vec<bool> {
        vec![
            true, true, true, true, true, true, true, true, false, false, false, false, false,
            false, false, false, true, true, true, true, false, false, false, false,
        ]
    }

    #[fixture]
    fn charge_control(
        simulation_time_iteration: SimulationTimeIteration,
        external_conditions: Arc<ExternalConditions>,
        external_sensor: ExternalSensor,
        charge_control_schedule: Vec<bool>,
    ) -> Arc<Control> {
        Arc::new(Control::Charge(
            ChargeControl::new(
                ControlLogicType::Automatic,
                charge_control_schedule,
                &simulation_time_iteration,
                0,
                1.,
                [1.0, 0.8].into_iter().map(Into::into).collect(),
                Some(22.),
                None,
                Some(external_conditions),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        ))
    }

    #[fixture]
    fn control() -> Arc<Control> {
        let mut schedule = vec![Some(21.), Some(21.), None, Some(21.)];
        schedule.extend(vec![None; 20]);

        Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            schedule,
            0,
            1.,
            Default::default(),
            Default::default(),
            1.,
        )))
    }

    fn create_elec_storage_heater(
        simulation_time: SimulationTime,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
        dry_core_min_output: Vec<[f64; 2]>,
        dry_core_max_output: Vec<[f64; 2]>,
        output_detailed_results: Option<bool>,
    ) -> Arc<ElecStorageHeater> {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "storage_heater").unwrap();

        let elec_storage_heater = ElecStorageHeater::new(
            3.5,
            2.5,
            10.0,
            ElectricStorageHeaterAirFlowType::FanAssisted,
            0.7,
            11.,
            1,
            21.,
            Arc::new(|| 20.),
            energy_supply_conn,
            &simulation_time.iter(),
            control,
            charge_control,
            dry_core_min_output,
            dry_core_max_output,
            external_conditions, // NOTE this is None in Python,
            0.,
            output_detailed_results,
        )
        .unwrap();

        elec_storage_heater.storage.write().set_state_of_charge(0.5);

        elec_storage_heater
    }

    #[fixture]
    fn elec_storage_heater(
        simulation_time: SimulationTime,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
    ) -> Arc<ElecStorageHeater> {
        create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            DRY_CORE_MAX_OUTPUT.to_vec(),
            None,
        )
    }

    #[rstest]
    fn test_initialisation(
        elec_storage_heater: Arc<ElecStorageHeater>,
        simulation_time_iteration: SimulationTimeIteration,
    ) {
        let energy = elec_storage_heater.demand_energy(1., &simulation_time_iteration);
        assert!(energy.unwrap() > 0.); // Should provide some energy

        // Test air flow type is set correctly by checking if fan energy is calculated
        assert_eq!(
            elec_storage_heater.air_flow_type,
            ElectricStorageHeaterAirFlowType::FanAssisted
        );
    }

    #[rstest]
    fn test_initialisation_invalid_soc_arrays(
        simulation_time: SimulationTime,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
    ) {
        let test_cases = [
            (
                vec![[0.0, 0.0], [0.5, 0.02], [0.3, 0.02], [1.0, 0.05]],
                vec![[0.0, 0.0], [0.5, 1.5], [0.7, 0.02], [1.0, 3.]],
                "shouldn't allow esh_min_output values in non-increasing order",
            ),
            (
                vec![[0.0, 0.0], [0.5, 0.02], [0.7, 0.02], [1.0, 0.05]],
                vec![[0.0, 0.0], [0.5, 1.5], [0.3, 0.02], [1.0, 3.]],
                "shouldn't allow esh_max_output values in non-increasing order",
            ),
            (
                vec![[0.0, 0.0], [0.5, 0.02], [0.7, 0.02], [0.9, 0.05]],
                vec![[0.0, 0.0], [0.5, 1.5], [0.7, 0.02], [1.0, 3.]],
                "shouldn't allow esh_min_output values not ending in 1.0",
            ),
            (
                vec![[0.0, 0.0], [0.5, 0.02], [0.7, 0.02], [1.0, 0.05]],
                vec![[0.0, 0.0], [0.5, 1.5], [0.7, 0.02], [0.9, 3.]],
                "shouldn't allow esh_max_output values not ending in 1.0",
            ),
            (
                vec![[0.2, 0.0], [0.5, 0.02], [0.7, 0.02], [1.0, 0.05]],
                vec![[0.0, 0.0], [0.5, 1.5], [0.7, 0.02], [1.0, 3.0]],
                "shouldn't allow esh_min_output values not starting at 1.0",
            ),
            (
                vec![[0.0, 0.0], [0.5, 0.02], [0.7, 0.02], [1.0, 0.05]],
                vec![[0.2, 0.0], [0.5, 1.5], [0.7, 0.02], [1.0, 3.0]],
                "shouldn't allow esh_max_output values not starting at 1.0",
            ),
            (
                vec![[0.0, 0.0], [1.0, 1.0]],
                vec![[0.0, 0.0], [1.0, 0.5]],
                "shouldn't allow any power_max values below power_min",
            ),
            (
                vec![[0.0, 0.0], [1.0, 1.0]],
                vec![[0.0, 0.0], [0.5, 0.4], [1.0, 1.0]],
                "shouldn't allow any power_max values below power_min",
            ),
        ];

        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "storage_heater").unwrap();

        for (esh_min_output, esh_max_output, debug_msg) in test_cases.iter() {
            assert!(
                ElecStorageHeater::new(
                    1.,
                    1.,
                    10.,
                    ElectricStorageHeaterAirFlowType::FanAssisted,
                    1.,
                    10.,
                    1,
                    21.,
                    Arc::new(|| 20.),
                    energy_supply_conn.clone(),
                    &simulation_time.iter(),
                    control.clone(),
                    charge_control.clone(),
                    esh_min_output.clone(),
                    esh_max_output.clone(),
                    external_conditions.clone(), // NOTE this is None in Python,
                    0.,
                    None,
                )
                .is_err(),
                "{}",
                debug_msg
            );
        }
    }

    #[rstest]
    fn test_initialisation_detailed_results(
        simulation_time: SimulationTime,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
    ) {
        let esh_max_output = vec![[0.0, 0.0], [0.5, 30.0], [1.0, 50.0]];
        let heater_with_detailed_results = create_elec_storage_heater(
            simulation_time,
            charge_control.clone(),
            control.clone(),
            external_conditions.clone(),
            DRY_CORE_MIN_OUTPUT.to_vec().clone(),
            esh_max_output.clone(),
            Some(true),
        );
        let heater_without_detailed_results = create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            esh_max_output,
            None,
        );

        assert!(heater_with_detailed_results.output_esh_results().is_some());
        assert!(heater_without_detailed_results
            .output_esh_results()
            .is_none());
    }

    #[rstest]
    fn test_temp_setpnt(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        let mut expected_setpoints = vec![Some(21.), Some(21.), None, Some(21.)];
        expected_setpoints.extend(vec![None; 20]);

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            assert_eq!(
                elec_storage_heater.temp_setpnt(&t_it),
                expected_setpoints[t_idx],
            );
        }
    }

    #[rstest]
    fn test_in_required_period(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        let mut expected_values = vec![true, true, false, true];
        expected_values.extend(vec![false; 20]);

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            assert_eq!(
                elec_storage_heater.in_required_period(&t_it),
                Some(expected_values[t_idx]),
            );
        }
    }

    #[rstest]
    fn test_frac_convective(elec_storage_heater: Arc<ElecStorageHeater>) {
        assert_eq!(elec_storage_heater.frac_convective(), 0.7);
    }

    #[rstest]
    #[ignore = "known issue"]
    fn test_energy_output_min(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        // Test minimum energy output calculation across all timesteps.

        let expected_min_energy_output = [
            0.019999999999999997,
            0.030419151282454364,
            0.0406826518009459,
            0.046014105108884346,
            0.04674323706081289,
            0.045800000000000014,
            0.046400000000000004,
            0.038,
            0.046886897763894056,
            0.03215061095009046,
            0.021233700713726503,
            0.01542628227371725,
            0.011428072222071864,
            0.008466125322130943,
            0.006271860545638567,
            0.004646308882570838,
            0.010432746398675731,
            0.01897277720547578,
            0.019999999999999997,
            0.019999999999999997,
            0.020056174317410178,
            0.014844117518306485,
            0.010996794239857175,
            0.008146626650951819,
        ]; // Actual minimum energy output for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let min_energy_output = elec_storage_heater.energy_output_min(&t_it).unwrap();
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);

            assert_relative_eq!(
                min_energy_output,
                expected_min_energy_output[t_idx],
                max_relative = 0.1
            );
        }
    }

    #[rstest]
    #[ignore = "known issue"]
    fn test_energy_output_max(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        // Test maximum energy output calculation across all timesteps.
        let expected_max_energy_output = [
            1.5,
            1.772121660521405,
            2.2199562136927717,
            2.5517202117781994,
            2.7913851590672585,
            2.7899999999999996,
            2.8200000000000003,
            2.4000000000000004,
            2.463423313846487,
            1.8249489529640162,
            1.3519554011630448,
            1.0015529506734968,
            0.7419686857505708,
            0.5496640579327374,
            0.40720123344887615,
            0.30166213975932143,
            0.6996897293886958,
            1.3814284569589004,
            1.5,
            1.5,
            1.3009346098448467,
            0.9637557636923015,
            0.713967931076402,
            0.5289205810615784,
        ]; // Expected max energy output for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let max_energy_output = elec_storage_heater.energy_output_max(&t_it).unwrap();
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let (energy, _, _, _) = max_energy_output;

            assert_relative_eq!(
                energy,
                expected_max_energy_output[t_idx],
                max_relative = 1e-1
            );
        }
    }

    #[rstest]
    fn test_energy_output_max_with_zero_event_single(
        simulation_time: SimulationTime,
        simulation_time_iterator: SimulationTimeIterator,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
    ) {
        let esh_max_output = vec![[0.0, 0.0], [0.5, 30.0], [1.0, 50.0]];
        let elec_storage_heater = create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            esh_max_output,
            None,
        );

        let expected_max_energy_output = 7.905696339321716;

        let expected_time_used = 1.0;

        let max_energy_output = elec_storage_heater
            .energy_output_max(&simulation_time_iterator.current_iteration())
            .unwrap();
        let _ =
            elec_storage_heater.demand_energy(5.0, &simulation_time_iterator.current_iteration());
        let (energy, time_used, _, _) = max_energy_output;

        assert_relative_eq!(energy, expected_max_energy_output, max_relative = 1e-2);

        assert_relative_eq!(
            time_used,
            expected_time_used,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[rstest]
    #[ignore = "known issue"]
    fn test_energy_output_max_with_zero_event(
        simulation_time: SimulationTime,
        simulation_time_iterator: SimulationTimeIterator,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
    ) {
        let elec_storage_heater = create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            DRY_CORE_MAX_OUTPUT.to_vec(),
            None,
        );

        // Test maximum energy output calculation across all timesteps.
        let expected_max_energy_output = [
            7.905696339321716,
            6.409423918606012,
            4.913144554873739,
            3.503509067065079,
            3.4999662250517263,
            3.5000272321002064,
            3.4999664471092706,
            3.5000359941220256,
            0.5819017195075272,
            0.0014406990152162622,
            7.97047184775446e-06,
            9.068336550842638e-08,
            0.0,
            0.0,
            0.0,
            0.0,
            2.9181276338320536,
            3.4985510038277714,
            3.5000309931442937,
            3.4999785113707134,
            0.5818631932808774,
            0.0014406043531798149,
            7.969521833559643e-06,
            9.066927928048405e-08,
        ];

        let expected_time_used = [
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.37318890798358917,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.37318890798725507,
        ];

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let max_energy_output = elec_storage_heater.energy_output_max(&t_it).unwrap();
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let (energy, time_used, _, _) = max_energy_output;

            assert_relative_eq!(
                energy,
                expected_max_energy_output[t_idx],
                max_relative = 1e-2
            );

            assert_relative_eq!(
                time_used,
                expected_time_used[t_idx],
                max_relative = EIGHT_DECIMAL_PLACES
            );
        }
    }

    #[rstest]
    fn test_electric_charge_automatic(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        // Test electric charge calculation across all timesteps.
        let expected_target_elec_charge = [
            0.5,
            1.0,
            0.99,
            0.98,
            0.95,
            0.93,
            0.9400000000000001,
            0.8,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.5,
            0.5,
            0.5,
            0.5,
            0.0,
            0.0,
            0.0,
            0.0,
        ]; // Expected target charge for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let target_elec_charge = elec_storage_heater.target_electric_charge(t_it).unwrap();
            assert_relative_eq!(target_elec_charge, expected_target_elec_charge[t_idx]);
        }
    }

    #[rstest]
    fn test_electric_charge_manual(
        external_conditions: Arc<ExternalConditions>,
        external_sensor: ExternalSensor,
        simulation_time: SimulationTime,
        control: Arc<Control>,
        charge_control_schedule: Vec<bool>,
    ) {
        let charge_control = Arc::new(Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                charge_control_schedule,
                &simulation_time.iter().current_iteration(),
                0,
                1.,
                [1.0, 0.8].into_iter().map(Into::into).collect(),
                Some(22.),
                None,
                Some(external_conditions.clone()),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        ));
        let heater = create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            DRY_CORE_MAX_OUTPUT.to_vec(),
            None,
        );
        let expected_target_elec_charge = [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let target_elec_charge = heater.target_electric_charge(t_it).unwrap();

            assert_relative_eq!(target_elec_charge, expected_target_elec_charge[t_idx]);
        }
    }

    #[rstest]
    fn test_electric_charge_celect(
        external_conditions: Arc<ExternalConditions>,
        external_sensor: ExternalSensor,
        simulation_time: SimulationTime,
        control: Arc<Control>,
        charge_control_schedule: Vec<bool>,
    ) {
        let charge_control = Arc::new(Control::Charge(
            ChargeControl::new(
                ControlLogicType::Celect,
                charge_control_schedule,
                &simulation_time.iter().current_iteration(),
                0,
                1.,
                [1.0, 0.8].into_iter().map(Into::into).collect(),
                Some(22.),
                None,
                Some(external_conditions.clone()),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        ));
        let heater = create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            DRY_CORE_MAX_OUTPUT.to_vec(),
            None,
        );
        let expected_target_elec_charge = [
            0.5, 1.0, 0.99, 0.98, 0.95, 0.93, 0.94, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let target_elec_charge = heater.target_electric_charge(t_it).unwrap();

            assert_relative_eq!(target_elec_charge, expected_target_elec_charge[t_idx]);
        }
    }

    #[rstest]
    fn test_electric_charge_hhrsh(
        external_conditions: Arc<ExternalConditions>,
        external_sensor: ExternalSensor,
        simulation_time: SimulationTime,
        control: Arc<Control>,
        charge_control_schedule: Vec<bool>,
    ) {
        let charge_control = Arc::new(Control::Charge(
            ChargeControl::new(
                ControlLogicType::Hhrsh,
                charge_control_schedule,
                &simulation_time.iter().current_iteration(),
                0,
                1.,
                [1.0, 0.8].into_iter().map(Into::into).collect(),
                Some(22.),
                None,
                Some(external_conditions.clone()),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        ));
        let heater = create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            DRY_CORE_MAX_OUTPUT.to_vec(),
            None,
        );
        let expected_target_elec_charge = [
            1., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 0., 0.,
            0., 0.,
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let target_elec_charge = heater.target_electric_charge(t_it).unwrap();

            assert_relative_eq!(target_elec_charge, expected_target_elec_charge[t_idx]);
        }
    }

    #[rstest]
    fn test_electric_charge_hhrsh_negative_heat_retention_ratio(
        external_conditions: Arc<ExternalConditions>,
        external_sensor: ExternalSensor,
        simulation_time: SimulationTime,
        control: Arc<Control>,
        charge_control_schedule: Vec<bool>,
    ) {
        let charge_control = Arc::new(Control::Charge(
            ChargeControl::new(
                ControlLogicType::Hhrsh,
                charge_control_schedule,
                &simulation_time.iter().current_iteration(),
                0,
                1.,
                [1.0, 0.8].into_iter().map(Into::into).collect(),
                Some(22.),
                None,
                Some(external_conditions.clone()),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        ));
        let heater = create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            DRY_CORE_MAX_OUTPUT.to_vec(),
            None,
        );

        heater.storage.write().set_heat_retention_ratio(-0.9);
        heater.storage.write().set_state_of_charge(0.5);

        let expected_target_elec_charge = [
            1., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 0., 0.,
            0., 0.,
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let target_elec_charge = heater.target_electric_charge(t_it).unwrap();

            assert_relative_eq!(target_elec_charge, expected_target_elec_charge[t_idx]);
        }
    }

    #[rstest]
    fn test_electric_charge_invalid_logic_type(
        external_conditions: Arc<ExternalConditions>,
        external_sensor: ExternalSensor,
        simulation_time: SimulationTime,
        control: Arc<Control>,
        charge_control_schedule: Vec<bool>,
    ) {
        let charge_control = Arc::new(Control::Charge(
            ChargeControl::new(
                ControlLogicType::HeatBattery,
                charge_control_schedule,
                &simulation_time.iter().current_iteration(),
                0,
                1.,
                [1.0, 0.8].into_iter().map(Into::into).collect(),
                Some(22.),
                None,
                Some(external_conditions.clone()),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        ));

        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "storage_heater").unwrap();

        let elec_storage_heater = ElecStorageHeater::new(
            3.5,
            2.5,
            10.0,
            ElectricStorageHeaterAirFlowType::FanAssisted,
            0.7,
            11.,
            1,
            21.,
            Arc::new(|| 20.),
            energy_supply_conn,
            &simulation_time.iter(),
            control,
            charge_control,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            DRY_CORE_MAX_OUTPUT.to_vec(),
            external_conditions, // NOTE this is None in Python,
            0.,
            None,
        );

        assert!(elec_storage_heater.is_err());
        assert_eq!(
            elec_storage_heater.unwrap_err().to_string(),
            "Control logic type HeatBattery is not valid for ElecStorageHeater."
        );
    }

    #[ignore = "fails probably because heat_retention_ratio cannot be set to None in Rust"]
    #[rstest]
    fn test_electric_charge_hhrsh_no_heat_retention_ratio(
        external_conditions: Arc<ExternalConditions>,
        external_sensor: ExternalSensor,
        simulation_time: SimulationTime,
        control: Arc<Control>,
        charge_control_schedule: Vec<bool>,
    ) {
        let charge_control = Arc::new(Control::Charge(
            ChargeControl::new(
                ControlLogicType::Hhrsh,
                charge_control_schedule,
                &simulation_time.iter().current_iteration(),
                0,
                1.,
                [1.0, 0.8].into_iter().map(Into::into).collect(),
                Some(22.),
                None,
                Some(external_conditions.clone()),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        ));

        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "storage_heater").unwrap();

        let heater = ElecStorageHeater::new(
            3.5,
            2.5,
            10.0,
            ElectricStorageHeaterAirFlowType::FanAssisted,
            0.7,
            11.,
            1,
            21.,
            Arc::new(|| 20.),
            energy_supply_conn,
            &simulation_time.iter(),
            control,
            charge_control,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            DRY_CORE_MAX_OUTPUT.to_vec(),
            external_conditions,
            0.,
            None,
        )
        .unwrap();

        heater.storage.write().set_heat_retention_ratio(0.);

        let result = heater
            .storage
            .read()
            .target_electric_charge(simulation_time.iter().current_iteration());

        assert!(result.is_err())
    }

    #[rstest]
    fn test_demand_energy(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        let expected_energy = [
            4.0,
            4.272121660521405,
            4.719956213692772,
            5.0,
            5.0,
            5.0,
            5.0,
            4.9,
            4.963423313846487,
            4.324948952964016,
            3.851955401163045,
            3.5015529506734966,
            3.241968685750571,
            3.0496640579327376,
            2.907201233448876,
            2.8016621397593213,
            3.199689729388696,
            3.8814284569589006,
            4.0,
            4.0,
            3.8009346098448464,
            3.4637557636923013,
            3.213967931076402,
            3.0289205810615782,
        ]; // Expected energy for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let energy_out = elec_storage_heater.demand_energy(5.0, &t_it).unwrap();
            assert_relative_eq!(energy_out, expected_energy[t_idx], max_relative = 1e-1);
        }
    }

    #[rstest]
    #[ignore = "known issue (energy_output)"]
    fn test_demand_energy_no_demand(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        let expected_energy = [
            0.019999999999999997,
            0.030419151282454364,
            0.04762390188197366,
            0.0488,
            0.04700000000000001,
            0.0458,
            0.046400000000000004,
            0.038,
            0.04932504653147142,
            0.04902998233007885,
            0.0487366832133454,
            0.0484451386224712,
            0.048155338061819486,
            0.04786727109853878,
            0.047580927362187296,
            0.047296296544359594,
            0.019999999999999997,
            0.019999999999999997,
            0.019999999999999997,
            0.019999999999999997,
            0.04701336839831551,
            0.046732132738611196,
            0.046452579440732555,
            0.04617469844073066,
        ]; // Expected energy for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let energy_out = elec_storage_heater.demand_energy(0., &t_it).unwrap();
            assert_relative_eq!(energy_out, expected_energy[t_idx], max_relative = 1e-1);
        }
    }

    #[rstest]
    fn test_demand_energy_detailed_results(
        simulation_time: SimulationTime,
        external_conditions: Arc<ExternalConditions>,
        control: Arc<Control>,
        charge_control: Arc<Control>,
    ) {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn = EnergySupply::connection(energy_supply, "storage_heater").unwrap();

        let heater = ElecStorageHeater::new(
            1.,
            1.,
            10.,
            ElectricStorageHeaterAirFlowType::FanAssisted,
            1.,
            10.,
            1,
            21.,
            Arc::new(|| 20.),
            energy_supply_conn,
            &simulation_time.iter(),
            control,
            charge_control,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            DRY_CORE_MAX_OUTPUT.to_vec(),
            external_conditions, // NOTE this is None in Python
            0.,
            Some(true),
        )
        .unwrap();

        let expected = [
            StorageHeaterDetailedResult {
                timestep_idx: 0,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.1360608282280595,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999999,
                energy_for_fan: 0.01,
                state_of_charge: 0.08639391717719404,
                final_soc: 0.08639391717719405,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 1,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.35997826917949083,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999998,
                energy_for_fan: 0.01,
                state_of_charge: 0.15039609025924494,
                final_soc: 0.15039609025924494,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 2,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.525860161130055,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999998,
                energy_for_fan: 0.01,
                state_of_charge: 0.1978100741462394,
                final_soc: 0.19781007414623947,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 3,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.6487484969563722,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999997,
                energy_for_fan: 0.01,
                state_of_charge: 0.23293522445060216,
                final_soc: 0.2329352244506022,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 4,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.7397864472877566,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999998,
                energy_for_fan: 0.01,
                state_of_charge: 0.2589565797218265,
                final_soc: 0.2589565797218265,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 5,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.8072290379637572,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999998,
                energy_for_fan: 0.01,
                state_of_charge: 0.2782336759254508,
                final_soc: 0.2782336759254508,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 6,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.8571917472765138,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999999,
                energy_for_fan: 0.01,
                state_of_charge: 0.29251450119779937,
                final_soc: 0.29251450119779937,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 7,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.894205037502244,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999998,
                energy_for_fan: 0.01,
                state_of_charge: 0.3030939974475749,
                final_soc: 0.303093997447575,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 8,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.7855640835379172,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.2245375890937832,
                final_soc: 0.2245375890937832,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 9,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.5819603347886747,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.16634155561491576,
                final_soc: 0.16634155561491573,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 10,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.4311269125454673,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.12322886436036903,
                final_soc: 0.12322886436036903,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 11,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.3193867248864036,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.09129019187172868,
                final_soc: 0.09129019187172868,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 12,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.23660753075068583,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.0676294387966601,
                final_soc: 0.06762943879666009,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 13,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.17528317748916877,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.050101121047743225,
                final_soc: 0.050101121047743225,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 14,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.1298529738329292,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.037115823664450306,
                final_soc: 0.037115823664450306,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 15,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.09619745025800532,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.027496078638649776,
                final_soc: 0.02749607863864978,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 16,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.20732566160442942,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999998,
                energy_for_fan: 0.01,
                state_of_charge: 0.10676351247820681,
                final_soc: 0.10676351247820684,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 17,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.41277253407174264,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999998,
                energy_for_fan: 0.01,
                state_of_charge: 0.16548625907103254,
                final_soc: 0.16548625907103257,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 18,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.5649711048613184,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999998,
                energy_for_fan: 0.01,
                state_of_charge: 0.20898914858490067,
                final_soc: 0.2089891485849007,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 19,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.6777226070926419,
                energy_instant: 1.0,
                energy_charged: 0.9999999999999998,
                energy_for_fan: 0.01,
                state_of_charge: 0.24121688787563644,
                final_soc: 0.24121688787563647,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 20,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.625190008414656,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.17869788703417083,
                final_soc: 0.17869788703417083,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 21,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.46315225418860434,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.1323826616153104,
                final_soc: 0.13238266161531043,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 22,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.34311168983701723,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.09807149263160868,
                final_soc: 0.09807149263160866,
                time_used_max: 1.0,
            },
            StorageHeaterDetailedResult {
                timestep_idx: 23,
                n_units: 1,
                energy_demand: 5.0,
                energy_delivered: 0.25418342247140796,
                energy_instant: 1.0,
                energy_charged: 0.0,
                energy_for_fan: 0.01,
                state_of_charge: 0.07265315038446787,
                final_soc: 0.07265315038446787,
                time_used_max: 1.0,
            },
        ];

        for t_it in simulation_time.iter() {
            heater.demand_energy(5., &t_it).unwrap();
        }

        let actual = heater.output_esh_results().unwrap();

        let assert_f64 = |expected: f64, actual: f64| {
            assert_relative_eq!(expected, actual, max_relative = 1e-6);
        };

        for t_idx in 0..24 {
            assert_eq!(expected[t_idx].timestep_idx, actual[t_idx].timestep_idx);
            assert_eq!(expected[t_idx].n_units, actual[t_idx].n_units);
            assert_f64(expected[t_idx].energy_demand, actual[t_idx].energy_demand);
            assert_f64(
                expected[t_idx].energy_delivered,
                actual[t_idx].energy_delivered,
            );
            assert_f64(expected[t_idx].energy_instant, actual[t_idx].energy_instant);
            assert_f64(expected[t_idx].energy_charged, actual[t_idx].energy_charged);
            assert_f64(expected[t_idx].energy_for_fan, actual[t_idx].energy_for_fan);
            assert_f64(
                expected[t_idx].state_of_charge,
                actual[t_idx].state_of_charge,
            );
            assert_f64(expected[t_idx].final_soc, actual[t_idx].final_soc);
            assert_f64(expected[t_idx].time_used_max, actual[t_idx].time_used_max);
        }
    }

    #[rstest]
    #[ignore = "known issue"]
    fn test_energy_for_fan(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        let expected_energy_for_fan = [
            0.003666666666666663,
            0.0034561513858793556,
            0.003133016027878955,
            0.0029047621124394848,
            0.0027330949708021268,
            0.0025983637642646956,
            0.0024893028132490476,
            0.002398927856611201,
            0.00243802553077851,
            0.0026117536650591346,
            0.002812153391189732,
            0.003045881325103026,
            0.0033220118367143004,
            0.003653245351058463,
            0.004057922556940931,
            0.004563550936407075,
            0.004342298713919709,
            0.0037302627836201,
            0.0036666666666666644,
            0.0036666666666666644,
            0.003861966010382251,
            0.004317143115057025,
            0.004894130588434899,
            0.005649492532217881,
        ]; // Expected energy for fan for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(0.5, &t_it);
            let energy_for_fan = elec_storage_heater.energy_for_fan();
            assert_relative_eq!(
                energy_for_fan,
                expected_energy_for_fan[t_idx],
                max_relative = EIGHT_DECIMAL_PLACES
            );
        }
    }

    #[rstest]
    #[ignore = "known issue"]
    fn test_energy_instant(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        let expected_energy_instant = [
            2.5,
            2.5,
            2.5,
            2.4482797882218006,
            2.208614840932741,
            2.2100000000000004,
            2.18,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
        ]; // Expected backup energy instant for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_instant = elec_storage_heater.energy_instant();
            assert_relative_eq!(
                energy_instant,
                expected_energy_instant[t_idx],
                max_relative = 0.1
            );
        }
    }

    #[rstest]
    #[ignore = "known issue"]
    fn test_energy_charged(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        let expected_energy_charged = [
            1.5,
            3.500000000000001,
            3.5,
            3.5,
            3.3397997646920756,
            2.7899999999999996,
            2.82,
            2.4000000000000004,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            3.500000000000001,
            2.738269398472182,
            1.5,
            1.5,
            0.0,
            0.0,
            0.0,
            0.0,
        ]; // Expected energy charged for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_charged = elec_storage_heater.energy_charged();
            assert_relative_eq!(
                energy_charged,
                expected_energy_charged[t_idx],
                max_relative = 0.1
            );
        }
    }

    #[rstest]
    #[ignore = "known issue"]
    fn test_energy_stored_delivered(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: Arc<ElecStorageHeater>,
    ) {
        let expected_energy_delivered = [
            1.5,
            1.772121660521405,
            2.219956213692772,
            2.5517202117781994,
            2.791385159067259,
            2.7899999999999996,
            2.82,
            2.4000000000000004,
            2.4634233138464876,
            1.8249489529640166,
            1.3519554011630448,
            1.0015529506734968,
            0.7419686857505706,
            0.5496640579327375,
            0.4072012334488761,
            0.30166213975932155,
            0.6996897293886956,
            1.3814284569589008,
            1.5,
            1.5,
            1.300934609844847,
            0.9637557636923016,
            0.713967931076402,
            0.5289205810615786,
        ]; // Expected energy stored delivered for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_delivered = elec_storage_heater.energy_delivered();
            assert_relative_eq!(
                energy_delivered,
                expected_energy_delivered[t_idx],
                max_relative = 1e-2
            );
        }
    }

    // skip test_invalid_air_flow_type as not replicable, Rust ensures only valid air flow types can be used

    #[rstest]
    fn test_damper_only(
        simulation_time: SimulationTime,
        external_conditions: Arc<ExternalConditions>,
        control: Arc<Control>,
        charge_control: Arc<Control>,
    ) {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn = EnergySupply::connection(energy_supply, "storage_heater").unwrap();

        let heater = ElecStorageHeater::new(
            3.5,
            2.5,
            10.0,
            ElectricStorageHeaterAirFlowType::DamperOnly,
            0.7,
            11.,
            1,
            21.,
            Arc::new(|| 20.),
            energy_supply_conn,
            &simulation_time.iter(),
            control,
            charge_control,
            DRY_CORE_MIN_OUTPUT.to_vec(),
            DRY_CORE_MAX_OUTPUT.to_vec(),
            external_conditions,
            0.,
            None,
        )
        .unwrap();

        heater.storage.write().set_state_of_charge(0.5);
        let _ = heater.demand_energy(1., &simulation_time.iter().current_iteration());

        assert_eq!(heater.energy_for_fan(), 0.)
    }

    #[rstest]
    fn test_protected_accessor_methods(elec_storage_heater: Arc<ElecStorageHeater>) {
        // skipped methods we did not port to Rust

        // Test getter methods
        assert_eq!(elec_storage_heater.storage.read().state_of_charge(), 0.5);
        assert_eq!(elec_storage_heater.storage.read().storage_capacity(), 10.);
        assert_eq!(elec_storage_heater.storage.read().n_units(), 1);
        assert_eq!(elec_storage_heater.storage.read().demand_met(), 0.);
        assert_eq!(elec_storage_heater.storage.read().demand_unmet(), 0.);

        // Test setter methods
        elec_storage_heater.storage.write().set_state_of_charge(0.8);
        assert_eq!(elec_storage_heater.storage.read().state_of_charge(), 0.8);

        elec_storage_heater.storage.write().set_demand_met(1.5);
        assert_eq!(elec_storage_heater.storage.read().demand_met(), 1.5);

        elec_storage_heater.storage.write().set_demand_unmet(0.5);
        assert_eq!(elec_storage_heater.storage.read().demand_unmet(), 0.5);

        // Test SOC clipping
        elec_storage_heater.storage.write().set_state_of_charge(1.5); // Above 1.0
        assert_eq!(elec_storage_heater.storage.read().state_of_charge(), 1.);

        elec_storage_heater.storage.write().set_state_of_charge(0.); // Below 1.0
        assert_eq!(elec_storage_heater.storage.read().state_of_charge(), 0.);
    }

    // TODO: test_invalid_charge_control_logic_error
    // TODO: test_elec_storage_heater_no_instant_power
    // TODO: test_elec_storage_energy_output_modes

    #[rstest]
    fn test_get_temp_for_charge_control(elec_storage_heater: Arc<ElecStorageHeater>) {
        let temp = elec_storage_heater.get_temp_for_charge_control().unwrap();
        assert_eq!(temp, 20.0);
    }

    #[rstest]
    fn test_get_zone_setpoint(elec_storage_heater: Arc<ElecStorageHeater>) {
        let setpoint = elec_storage_heater.get_zone_setpoint();
        assert_eq!(setpoint, 21.);
    }

    // TODO: test_output_mode_enum_values
    // TODO: test_heat_storage_dry_core_boundary_soc_values
    // TODO: test_heat_storage_dry_core_zero_capacity
}
