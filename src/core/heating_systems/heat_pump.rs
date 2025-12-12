/// This module provides objects to represent heat pumps and heat pump test data.
/// The calculations are based on the DAHPSE method developed for generating PCDB
/// entries for SAP 2012 and SAP 10. DAHPSE was based on a draft of
/// BS EN 15316-4-2:2017 and is described in the SAP calculation method CALCM-01.
use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::boiler::{Boiler, BoilerServiceWaterCombi};
use crate::core::heating_systems::boiler::{BoilerServiceSpace, BoilerServiceWaterRegular};
use crate::core::heating_systems::common::HeatingServiceType;
use crate::core::material_properties::{MaterialProperties, WATER};
use crate::core::schedule::{expand_numeric_schedule, reject_nulls};
use crate::core::units::{
    celsius_to_kelvin, kelvin_to_celsius, BelowAbsoluteZeroError, HOURS_PER_DAY,
    KILOJOULES_PER_KILOWATT_HOUR, SECONDS_PER_MINUTE, WATTS_PER_KILOWATT,
};
use crate::corpus::{ResultParamValue, ResultsAnnual, ResultsPerTimestep, TempInternalAirFn};
use crate::external_conditions::ExternalConditions;
use crate::input::{
    BoilerCostScheduleHybrid, HeatPumpBackupControlType, HeatPumpHotWaterOnlyTestDatum,
    HeatPumpHotWaterTestData, HeatPumpSinkType, HeatPumpSourceType,
    HeatPumpTestDatum as HeatPumpTestDatumInput, HeatSourceWetDetails, HotWaterSourceDetails,
};
use crate::simulation_time::SimulationTimeIteration;
use crate::statistics::np_interp;
use anyhow::{anyhow, bail};
use derivative::Derivative;
use indexmap::IndexMap;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use parking_lot::{Mutex, RwLock};
use polyfit_rs::polyfit_rs::polyfit;
use serde::{Deserialize, Serialize};
use smartstring::alias::String;
use std::collections::HashMap;
use std::fmt::{Debug, Formatter};
use std::sync::Arc;


const N_EXER: f64 = 3.0;

impl HeatPumpSourceType {
    pub fn is_exhaust_air(&self) -> bool {
        matches!(
            self,
            HeatPumpSourceType::ExhaustAirMEV
                | HeatPumpSourceType::ExhaustAirMVHR
                | HeatPumpSourceType::ExhaustAirMixed
        )
    }

    pub fn source_fluid_is_air(&self) -> bool {
        matches!(
            self,
            HeatPumpSourceType::OutsideAir
                | HeatPumpSourceType::ExhaustAirMEV
                | HeatPumpSourceType::ExhaustAirMVHR
                | HeatPumpSourceType::ExhaustAirMixed
        )
    }

    pub fn source_fluid_is_water(&self) -> bool {
        matches!(
            self,
            HeatPumpSourceType::Ground
                | HeatPumpSourceType::WaterGround
                | HeatPumpSourceType::WaterSurface
                | HeatPumpSourceType::HeatNetwork
        )
    }
}

#[derive(Clone, Debug)]
pub(crate) enum HeatPumpEmitterType {
    RadiatorsUfh,
    FanCoils,
    WarmAir,
}

/// Calculate Carnot CoP based on source and outlet temperatures (in Kelvin)
fn carnot_cop(temp_source: f64, temp_outlet: f64, temp_diff_limit_low: Option<f64>) -> f64 {
    let mut temp_diff = temp_outlet - temp_source;
    if let Some(low_limit) = temp_diff_limit_low {
        temp_diff = max_of_2(temp_diff, low_limit);
    }
    temp_outlet / temp_diff
}

/// Interpolate between test data records for different air flow rates
///
/// Arguments:
/// * `throughput_exhaust_air` - throughput (m3 / h) of exhaust air
/// * `test_data`
///        - list of records of heat pump test data, each with the following elements:
///                - air_flow_rate
///                - test_letter
///                - capacity
///                - cop
///                - design_flow_temp (in Celsius)
///                - temp_outlet (in Celsius)
///                - temp_source (in Celsius)
///                - temp_test (in Celsius)
fn interpolate_exhaust_air_heat_pump_test_data(
    throughput_exhaust_air: f64,
    test_data: &Vec<HeatPumpTestDatumInput>,
    source_type: HeatPumpSourceType,
) -> anyhow::Result<(f64, Vec<HeatPumpTestDatum>)> {
    // split test records into different lists by air flow rate
    let mut test_data_by_air_flow_rate: HashMap<OrderedFloat<f64>, Vec<HeatPumpTestDatum>> =
        Default::default();
    for test_data_record in test_data {
        if test_data_record.air_flow_rate.is_none() {
            continue;
        }
        test_data_by_air_flow_rate
            .entry(OrderedFloat(test_data_record.air_flow_rate.unwrap()))
            .or_default()
            .push(HeatPumpTestDatum::from(test_data_record.clone()));
    }

    // check that all lists have same combination of design flow temp and test letter
    let mut fixed_temps_and_test_letters: Option<Vec<TestDatumTempsAndTestLetters>> = None;
    for (air_flow_rate, test_data_record_list) in &test_data_by_air_flow_rate {
        // find and save all the combinations of design flow temp and test letter for this air flow rate
        let mut air_flowrate_flow_temps_test_letters: Vec<TestDatumTempsAndTestLetters> =
            Default::default();
        for test_data_record in test_data_record_list {
            air_flowrate_flow_temps_test_letters.push(TestDatumTempsAndTestLetters {
                air_flow_rate: air_flow_rate.0,
                design_flow_temp: test_data_record.design_flow_temp,
                test_letter: test_data_record.test_letter,
                temp_outlet: test_data_record.temp_outlet,
                temp_source: test_data_record.temp_source,
                temp_test: test_data_record.temp_test,
            })
        }

        match fixed_temps_and_test_letters {
            None => {
                fixed_temps_and_test_letters = Some(air_flowrate_flow_temps_test_letters);
            }
            Some(ref fixed_temps_and_test_letters) => {
                if fixed_temps_and_test_letters != &air_flowrate_flow_temps_test_letters {
                    bail!("In heat pump test data, fixed temps and test_letters were not consistent for one air_flow_temp value");
                }
            }
        }
    }

    // Construct test data records interpolated by air flow rate
    let air_flow_rates_ordered = &test_data_by_air_flow_rate
        .keys()
        .sorted()
        .map(|f| f.0)
        .collect::<Vec<f64>>();
    let test_data_interp_by_air_flow_rate = fixed_temps_and_test_letters
        .ok_or(anyhow!("Non-empty test data was expected"))?
        .iter()
        .map(|datum| {
            // create lists of test data values ordered by air flow rate
            let mut capacity_list: Vec<f64> = Default::default();
            let mut cop_list: Vec<f64> = Default::default();
            let mut ext_air_ratio_list: Vec<f64> = Default::default();
            for air_flow_rate in air_flow_rates_ordered {
                for test_record in &test_data_by_air_flow_rate[&OrderedFloat(*air_flow_rate)] {
                    if test_record.design_flow_temp == datum.design_flow_temp
                        && test_record.test_letter == datum.test_letter
                    {
                        capacity_list.push(test_record.capacity);
                        cop_list.push(test_record.cop);
                        if source_type == HeatPumpSourceType::ExhaustAirMixed {
                            ext_air_ratio_list.push(test_record.eahp_mixed_ext_air_ratio.ok_or(anyhow!("An eahp_mixed_ext_air_ratio value was expected for a heat pump with source type ExhaustAirMixed."))?);
                        }
                    }
                }
            }

            let capacity = np_interp(
                throughput_exhaust_air, air_flow_rates_ordered,
                &capacity_list,
            );
            let cop = np_interp(throughput_exhaust_air, air_flow_rates_ordered, &cop_list);
            let ext_air_ratio =
                (source_type == HeatPumpSourceType::ExhaustAirMixed).then(|| np_interp(
                    throughput_exhaust_air, air_flow_rates_ordered,
                    &ext_air_ratio_list,
                ));

            Ok(HeatPumpTestDatum {
                air_flow_rate: Some(datum.air_flow_rate),
                test_letter: datum.test_letter,
                capacity,
                cop,
                design_flow_temp: datum.design_flow_temp,
                temp_outlet: datum.temp_outlet,
                temp_source: datum.temp_source,
                temp_test: datum.temp_test,
                eahp_mixed_ext_air_ratio: ext_air_ratio,
            })
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    let lowest_air_flow_rate_in_test_data = &air_flow_rates_ordered
        .iter()
        .max_by(|a, b| a.total_cmp(b).reverse())
        .ok_or(anyhow!("Non-empty test data was expected"))?;

    Ok((
        lowest_air_flow_rate_in_test_data.to_owned().to_owned(),
        test_data_interp_by_air_flow_rate,
    ))
}

#[derive(Clone, Debug)]
pub(crate) struct HeatPumpTestDatum {
    air_flow_rate: Option<f64>,
    test_letter: TestLetter,
    capacity: f64,
    cop: f64,
    design_flow_temp: f64,
    temp_outlet: f64,
    temp_source: f64,
    temp_test: f64,
    eahp_mixed_ext_air_ratio: Option<f64>,
}

impl From<HeatPumpTestDatumInput> for HeatPumpTestDatum {
    fn from(value: HeatPumpTestDatumInput) -> Self {
        let HeatPumpTestDatumInput {
            air_flow_rate,
            test_letter,
            capacity,
            cop,
            design_flow_temp,
            temp_outlet,
            temp_source,
            temp_test,
            eahp_mixed_ext_air_ratio,
        } = value;
        Self {
            air_flow_rate,
            test_letter,
            capacity,
            cop,
            design_flow_temp,
            temp_outlet,
            temp_source,
            temp_test,
            eahp_mixed_ext_air_ratio,
        }
    }
}

//part of the thermal losses transmitted to the room
//Taken from Storage tank module: Method A from BS EN 15316-5:2017
const F_STO_M: f64 = 0.75;

// Reference conditions for buffer tank thermal loss calculations
const TEMP_SET_REF: f64 = 65.; // Reference temperature setpnt in °C
const TEMP_AMB_REF: f64 = 20.; // Reference ambient temperature in °C

// Temperature split factor for flow/return calculations
const TEMP_SPLIT_FACTOR_FLOW: f64 = 0.5;

#[derive(Clone, Copy, Debug, Default)]
pub struct BufferTankEmittersData {
    pub temp_emitter_req: f64,
    pub power_req_from_buffer_tank: f64,
    pub design_flow_temp: f64,
    pub target_flow_temp: f64,
    pub temp_rm_prev: f64,
    pub variable_flow: bool,
    pub temp_diff_emit_design: f64,
    pub min_flow_rate: f64,
    pub max_flow_rate: f64,
}

impl BufferTankEmittersData {
    fn merge_result(self, result: &BufferTankServiceResult) -> BufferTankEmittersDataWithResult {
        BufferTankEmittersDataWithResult {
            data: self,
            result: result.clone(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct BufferTankEmittersDataWithResult {
    pub data: BufferTankEmittersData,
    pub result: BufferTankServiceResult,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BufferTankServiceResult {
    _service_name: String,
    _power_req_from_buffer_tank: f64,
    _temp_emitter_req: f64,
    _buffer_emitter_circ_flow_rate: f64,
    flow_temp_increase_due_to_buffer: f64,
    pump_power_at_flow_rate: f64,
    heat_loss_buffer_kwh: f64,
}

/// A type for Buffer Tanks linked to heat pumps.
#[derive(Clone, Debug)]
pub struct BufferTank {
    daily_losses: f64,
    volume: f64,
    pump_fixed_flow_rate: f64,
    pump_power_at_flow_rate: f64,
    number_of_zones: usize,
    simulation_timestep: f64,
    _contents: MaterialProperties,
    service_results: Vec<BufferTankServiceResult>,
    _buffer_emitter_circ_flow_rate: f64,
    _hp_buffer_circ_flow_rate: f64,
    _print_warning: bool,
    cp: f64,
    rho: f64,
    q_heat_loss_buffer_rbl: f64,
    track_buffer_loss: f64,
    temp_average_buffer: f64,
    detailed_results: Option<Vec<Vec<BufferTankServiceResult>>>,
}

impl BufferTank {
    pub fn new(
        daily_losses: f64,
        volume: f64,
        pump_fixed_flow_rate: f64,
        pump_power_at_flow_rate: f64,
        number_of_zones: usize,
        simulation_timestep: f64,
        initial_temp: f64,
        contents: MaterialProperties,
        output_detailed_results: Option<bool>,
    ) -> Self {
        let output_detailed_results = output_detailed_results.unwrap_or(false);
        let (cp, rho) = (contents.specific_heat_capacity_kwh(), contents.density());

        Self {
            daily_losses,
            volume,
            pump_fixed_flow_rate,
            pump_power_at_flow_rate,
            number_of_zones,
            simulation_timestep,
            _contents: contents,
            service_results: Default::default(),
            _buffer_emitter_circ_flow_rate: 0.0,
            _hp_buffer_circ_flow_rate: 0.0,
            _print_warning: true,
            cp,
            rho,
            q_heat_loss_buffer_rbl: 0.0,
            track_buffer_loss: 0.0,
            // Initialisation of buffer tank temperature needed in case no power required by emitters in first step
            temp_average_buffer: initial_temp,
            detailed_results: output_detailed_results.then(Default::default),
        }
    }

    pub fn update_buffer_loss(&mut self, buffer_loss: f64) {
        self.track_buffer_loss = buffer_loss;
    }

    pub fn buffer_loss(&self) -> f64 {
        self.track_buffer_loss
    }

    /// Thermal losses are calculated with respect to the impact of the temperature set point
    fn thermal_losses(&mut self, temp_average_buffer: f64, temp_rm_prev: f64) -> f64 {
        // Specific loss in W/K
        let h_buffer_ls = (WATTS_PER_KILOWATT as f64 * self.daily_losses)
            / (HOURS_PER_DAY as f64 * (TEMP_SET_REF - TEMP_AMB_REF));
        // Absolute loss in W
        let heat_loss_buffer_w = (temp_average_buffer - temp_rm_prev) * h_buffer_ls;
        let heat_loss_buffer_kwh = heat_loss_buffer_w / WATTS_PER_KILOWATT as f64
            * self.simulation_timestep
            / self.number_of_zones as f64;

        // TODO (from Python): update when zoning sorted out. Hopefully then there will be no need to divide by the number of zones.

        // recoverable heat losses from buffer - kWh
        self.q_heat_loss_buffer_rbl = heat_loss_buffer_kwh * F_STO_M;

        heat_loss_buffer_kwh
    }

    /// Return the recoverable heat losses as internal gain for the current timestep in W
    pub fn internal_gains(&self) -> f64 {
        // multiply by number of zones because this is only called once (not per zone like thermal_loss) for internal gains
        self.q_heat_loss_buffer_rbl * WATTS_PER_KILOWATT as f64 / self.simulation_timestep
            * self.number_of_zones as f64
    }

    pub fn calc_buffer_tank(
        &mut self,
        service_name: &str,
        emitters_data_for_buffer_tank: BufferTankEmittersData,
    ) -> anyhow::Result<&[BufferTankServiceResult]> {
        let temp_rm_prev = emitters_data_for_buffer_tank.temp_rm_prev;

        let result_service_name = String::from([service_name, "_buffer_tank"].concat());

        if emitters_data_for_buffer_tank.power_req_from_buffer_tank > 0.0 {
            let temp_emitter_req = emitters_data_for_buffer_tank.temp_emitter_req;
            self.temp_average_buffer = temp_emitter_req;

            // call to calculate thermal losses
            let heat_loss_buffer_kwh = self.thermal_losses(temp_emitter_req, temp_rm_prev);

            // We are calculating the delta T needed to give the heat output required from the emitters
            // E = m C dT rearranged to make delta T the subject
            let delta_t_buffer = (emitters_data_for_buffer_tank.power_req_from_buffer_tank
                / (self.pump_fixed_flow_rate / SECONDS_PER_MINUTE as f64))
                / (self.cp * self.rho * KILOJOULES_PER_KILOWATT_HOUR as f64);

            let buffer_flow_temp = temp_emitter_req + TEMP_SPLIT_FACTOR_FLOW * delta_t_buffer;
            let buffer_return_temp =
                temp_emitter_req - (1. - TEMP_SPLIT_FACTOR_FLOW) * delta_t_buffer;
            let theoretical_hp_return_temp = buffer_return_temp;

            let mut flag = true;
            let mut hp_flow;
            let mut theoretical_hp_flow_temp: Option<f64> = None;
            if emitters_data_for_buffer_tank.variable_flow {
                let delta_t_hp_to_buffer = emitters_data_for_buffer_tank.temp_diff_emit_design;
                theoretical_hp_flow_temp = Some(theoretical_hp_return_temp + delta_t_hp_to_buffer);
                hp_flow = (emitters_data_for_buffer_tank.power_req_from_buffer_tank
                    + heat_loss_buffer_kwh / self.simulation_timestep)
                    / (delta_t_hp_to_buffer
                        * self.cp
                        * self.rho
                        * KILOJOULES_PER_KILOWATT_HOUR as f64);
                let max_flow_rate = emitters_data_for_buffer_tank.max_flow_rate;
                let min_flow_rate = emitters_data_for_buffer_tank.min_flow_rate;

                if hp_flow > max_flow_rate {
                    hp_flow = max_flow_rate
                } else if hp_flow < min_flow_rate {
                    hp_flow = min_flow_rate
                } else {
                    flag = false;
                };
            } else {
                hp_flow = emitters_data_for_buffer_tank.min_flow_rate;
            }

            if flag {
                theoretical_hp_flow_temp = Some(
                    theoretical_hp_return_temp
                        + ((emitters_data_for_buffer_tank.power_req_from_buffer_tank
                            + heat_loss_buffer_kwh / self.simulation_timestep)
                            / hp_flow)
                            / (self.cp * self.rho * KILOJOULES_PER_KILOWATT_HOUR as f64),
                );
            }

            // We are currently assuming that, by design, buffer tanks always work with fix pumps on the
            // tank-emitter side with higher flow than the flow in the hp-tank side.
            if hp_flow > self.pump_fixed_flow_rate / SECONDS_PER_MINUTE as f64 {
                bail!("Heat pump buffer tank flow is greater than buffer tank emitter flow. Calculation aborted")
            }

            let flow_temp_increase_due_to_buffer = max_of_2(
                0.,
                theoretical_hp_flow_temp.expect("Heat pump calc_buffer_tank function expects theoretical_hp_flow_temp to be set")
                    - buffer_flow_temp,
            );

            // TODO (from Python) Consider circumstances where the flow rate on the hp-buffer loop exceed the flow rate in the
            //      buffer-emitter loop. We think in this circumstances the flow temperature would match rather than
            //      the return temperature, changing the logic of the methodology.
            //      This would happen in situation of very cold weather which exceed the design conditions for the
            //      buffer-emitters loop sizing/flow rate

            // If detailed results are to be output, save the results from the current timestep
            self.service_results.push(BufferTankServiceResult {
                _service_name: result_service_name,
                _power_req_from_buffer_tank: emitters_data_for_buffer_tank
                    .power_req_from_buffer_tank,
                _temp_emitter_req: temp_emitter_req,
                _buffer_emitter_circ_flow_rate: self.pump_fixed_flow_rate,
                flow_temp_increase_due_to_buffer,
                pump_power_at_flow_rate: self.pump_power_at_flow_rate,
                heat_loss_buffer_kwh,
            });

            if let Some(detailed_results) = self.detailed_results.as_mut() {
                detailed_results.push(self.service_results.clone());
            }
        } else {
            // call to calculate cool down losses
            let heat_loss_buffer_kwh = self.thermal_losses(self.temp_average_buffer, temp_rm_prev);
            let heat_capacity_buffer =
                self.volume * (self.rho * self.cp * KILOJOULES_PER_KILOWATT_HOUR as f64);
            let temp_loss =
                heat_loss_buffer_kwh / (heat_capacity_buffer / KILOJOULES_PER_KILOWATT_HOUR as f64);
            let new_temp_average_buffer = self.temp_average_buffer - temp_loss;

            self.temp_average_buffer = new_temp_average_buffer;

            self.service_results.push(BufferTankServiceResult {
                _service_name: result_service_name,
                _power_req_from_buffer_tank: 0.0,
                _temp_emitter_req: emitters_data_for_buffer_tank.temp_emitter_req,
                _buffer_emitter_circ_flow_rate: self.pump_fixed_flow_rate,
                flow_temp_increase_due_to_buffer: 0.0,
                pump_power_at_flow_rate: 0.0,
                heat_loss_buffer_kwh,
            });
            if let Some(detailed_results) = self.detailed_results.as_mut() {
                detailed_results.push(self.service_results.clone());
            }
        }

        Ok(&self.service_results)
    }

    #[cfg(test)]
    pub(crate) fn set_track_buffer_loss(&mut self, buffer_loss: f64) {
        self.track_buffer_loss = buffer_loss;
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum TestLetter {
    A,
    B,
    C,
    D,
    F,
}

impl std::fmt::Display for TestLetter {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                TestLetter::A => "A",
                TestLetter::B => "B",
                TestLetter::C => "C",
                TestLetter::D => "D",
                TestLetter::F => "F",
            }
        )
    }
}

impl TestLetter {
    fn is_non_bivalent(&self) -> bool {
        [TestLetter::A, TestLetter::B, TestLetter::C, TestLetter::D].contains(self)
    }
}

#[derive(Derivative)]
#[derivative(PartialEq, Debug)]
pub struct TestDatumTempsAndTestLetters {
    #[derivative(PartialEq = "ignore")]
    // we need to ignore air flow rate when performing comparisons
    pub air_flow_rate: f64,
    pub design_flow_temp: f64,
    pub(crate) test_letter: TestLetter,
    pub temp_outlet: f64,
    pub temp_source: f64,
    pub temp_test: f64,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CompleteHeatPumpTestDatum {
    pub air_flow_rate: Option<f64>,
    pub(crate) test_letter: TestLetter,
    pub capacity: f64,
    pub cop: f64,
    pub design_flow_temp: f64,
    pub temp_outlet: f64,
    pub temp_source: f64,
    pub temp_test: f64,
    pub carnot_cop: f64,
    pub exergetic_eff: f64,
    pub theoretical_load_ratio: f64,
}

impl CompleteHeatPumpTestDatum {
    pub fn data_item(&self, item: &DatumItem) -> f64 {
        match item {
            DatumItem::CarnotCop => self.carnot_cop,
            DatumItem::TempOutlet => self.temp_outlet,
            DatumItem::TempSource => self.temp_source,
            DatumItem::Capacity => self.capacity,
        }
    }
}

pub enum DatumItem {
    CarnotCop,
    TempOutlet,
    TempSource,
    Capacity,
}

impl HeatPumpTestDatum {
    pub fn complete(
        &self,
        carnot_cop: f64,
        exergetic_eff: f64,
        theoretical_load_ratio: f64,
    ) -> CompleteHeatPumpTestDatum {
        let HeatPumpTestDatum {
            air_flow_rate,
            test_letter,
            capacity,
            cop,
            design_flow_temp,
            temp_outlet,
            temp_source,
            temp_test,
            ..
        } = self;
        CompleteHeatPumpTestDatum {
            air_flow_rate: *air_flow_rate,
            test_letter: *test_letter,
            capacity: *capacity,
            cop: *cop,
            design_flow_temp: *design_flow_temp,
            temp_outlet: *temp_outlet,
            temp_source: *temp_source,
            temp_test: *temp_test,
            carnot_cop,
            exergetic_eff,
            theoretical_load_ratio,
        }
    }
}

/// An object to represent EN 14825 test data for a heat pump.
///
/// This object stores the data and provides functions to look up values from
/// the correct data records for the conditions being modelled.
///
/// NB. OrderedFloat values are used as keys for the test data as,
/// unlike f64, this is a representation of a float that is both Hash + Eq and so can be
/// used as a key in a HashMap.
#[derive(Clone, Debug)]
pub(crate) struct HeatPumpTestData {
    test_data: HashMap<OrderedFloat<f64>, Vec<CompleteHeatPumpTestDatum>>,
    design_flow_temps: Vec<OrderedFloat<f64>>,
    average_cap: Vec<f64>,
    temp_spread_test_conditions: Vec<f64>,
    regression_coeffs: HashMap<OrderedFloat<f64>, Vec<f64>>,
}

const TEST_LETTERS_NON_BIVALENT: [char; 4] = ['A', 'B', 'C', 'D'];
const TEST_LETTERS_ALL: [char; 5] = ['A', 'B', 'C', 'D', 'F'];

impl HeatPumpTestData {
    pub(crate) fn new(data: Vec<HeatPumpTestDatum>) -> anyhow::Result<Self> {
        // keyed by design flow temp
        let mut test_data: HashMap<OrderedFloat<f64>, Vec<HeatPumpTestDatum>> = Default::default();
        let mut design_flow_temps: Vec<OrderedFloat<f64>> = Default::default();

        // variable to count duplicate records for each design flow temp
        let mut dupl: HashMap<OrderedFloat<f64>, usize> = Default::default();

        for datum in data {
            let mut saved_datum = datum.clone();
            let design_flow_temp = OrderedFloat(datum.design_flow_temp);

            // When a new design flow temp is encountered, add it to lists/ maps
            if !dupl.contains_key(&design_flow_temp) {
                if !design_flow_temps.contains(&design_flow_temp) {
                    design_flow_temps.push(design_flow_temp);
                }
                test_data.entry(design_flow_temp).or_default();
            }

            let mut duplicate = false;
            for d in test_data.get(&design_flow_temp).unwrap() {
                if datum == *d {
                    duplicate = true;
                    // Increment count of number of duplicates for this design flow temp
                    // Handle records with same inlet temp
                    // Cannot process a row at the same inlet temperature (div
                    // by zero error during interpolation), so we add a tiny
                    // amount to the temperature (to 10DP) for each duplicate
                    // found.
                    saved_datum.temp_test += 0.0000000001;
                    saved_datum.temp_source += 0.0000000001;
                }
            }
            // This increment has to be after loop to avoid multiple-counting
            // when there are 3 or more duplicates. E.g. if there are already 2
            // records that are the same, then when adding a third that is the
            // same, we only want to increment the counter by 1 (for the record
            // we are adding) and not 2 (the number of existing records the new
            // record duplicates).
            if duplicate {
                *dupl.entry(design_flow_temp).or_default() += 1;
            }

            test_data
                .entry(design_flow_temp)
                .or_default()
                .push(saved_datum);
        }

        // Check the number of test records is as expected
        // - Up to 4 design flow temps (EN 14825:2018 defines these as 35, 45, 55, 65)
        // - 4 or 5 distinct records for each flow temp
        if design_flow_temps.is_empty() {
            bail!("No test data provided for heat pump performance");
        } else if design_flow_temps.len() > 4 {
            bail!("Warning: Test data for a maximum of 4 design flow temperatures is expected. {} have been provided.", design_flow_temps.len());
        }
        for (design_flow_temp, data) in test_data.iter() {
            if dupl.contains_key(design_flow_temp) {
                if (data.len() - dupl.get(design_flow_temp).unwrap()) != 4 {
                    bail!("Expected 4 distinct records for each design flow temperature");
                }
            } else if data.len() != 5 {
                bail!("Expected 5 records for each design flow temperature");
            }
        }

        // Check if test letters ABCDE are present as expected
        let mut test_letter_vec: Vec<char> = Default::default();
        for temperature in &design_flow_temps {
            for test_data in &test_data[temperature] {
                for test_letter in test_data.test_letter.to_string().chars() {
                    test_letter_vec.push(test_letter);
                }
                if test_letter_vec.len() == 5 {
                    for test_letter_check in TEST_LETTERS_ALL {
                        if !test_letter_vec.contains(&test_letter_check) {
                            bail!("Expected test letter {test_letter_check} in {temperature} degree temp data");
                        }
                    }
                    test_letter_vec = Default::default();
                }
            }
        }

        design_flow_temps.sort_by(|a, b| a.total_cmp(b));
        for (_, data) in test_data.iter_mut() {
            data.sort_by(|a, b| a.temp_test.total_cmp(&b.temp_test));
        }

        let average_cap = average_capacity(&design_flow_temps, &test_data);
        let temp_spread_test_conditions = init_temp_spread_test_conditions(&design_flow_temps)?;
        let regression_coeffs = init_regression_coeffs(&design_flow_temps, &test_data)?;

        // Calculate derived variables for each data record which are not time-dependent
        let test_data: HashMap<OrderedFloat<f64>, Vec<CompleteHeatPumpTestDatum>> = test_data
            .iter()
            .map(|(design_flow_temp, data)| {
                let (carnot_cops, exergetic_effs) =
                    data.iter().fold((vec![], vec![]), |mut acc, datum| {
                        // Get the source and outlet temperatures from the test record
                        let temp_source = celsius_to_kelvin(datum.temp_source)
                            .expect("Invalid temperature encountered");
                        let temp_outlet = celsius_to_kelvin(datum.temp_outlet)
                            .expect("Invalid temperature encountered");

                        // Calculate the Carnot CoP
                        let carnot_cop = carnot_cop(temp_source, temp_outlet, None);
                        acc.0.push(carnot_cop);
                        acc.1.push(datum.cop / carnot_cop);
                        acc
                    });

                let temp_source_cld = celsius_to_kelvin(data[0].temp_source)
                    .expect("Invalid temperature not expected");
                let temp_outlet_cld = celsius_to_kelvin(data[0].temp_outlet)
                    .expect("Invalid temperature not expected");
                let carnot_cop_cld = carnot_cops[0];

                let theoretical_load_ratios =
                    data.iter().enumerate().fold(vec![], |mut acc, (i, datum)| {
                        // Get the source and outlet temperatures from the test record
                        let temp_source = celsius_to_kelvin(datum.temp_source)
                            .expect("Invalid temperature not expected");
                        let temp_outlet = celsius_to_kelvin(datum.temp_outlet)
                            .expect("Invalid temperature not expected");

                        let theoretical_load_ratio = (carnot_cops[i] / carnot_cop_cld)
                            * (temp_outlet_cld * temp_source / (temp_source_cld * temp_outlet))
                                .powf(N_EXER);
                        acc.push(theoretical_load_ratio);
                        acc
                    });

                let complete_data: Vec<CompleteHeatPumpTestDatum> = data
                    .iter()
                    .enumerate()
                    .map(|(i, datum)| {
                        datum.complete(
                            carnot_cops[i],
                            exergetic_effs[i],
                            theoretical_load_ratios[i],
                        )
                    })
                    .collect();
                (*design_flow_temp, complete_data)
            })
            .collect();

        Ok(Self {
            test_data,
            design_flow_temps,
            average_cap,
            temp_spread_test_conditions,
            regression_coeffs,
        })
    }

    /// Return average capacity for tests A-D, interpolated between design flow temps
    pub fn average_capacity(
        &self,
        design_flow_temp_op_cond: f64,
    ) -> Result<f64, BelowAbsoluteZeroError> {
        if self.design_flow_temps.len() == 1 {
            return Ok(self.average_cap[0]);
        }

        let flow_temp = kelvin_to_celsius(design_flow_temp_op_cond)?;
        Ok(np_interp(
            flow_temp,
            &self
                .design_flow_temps
                .iter()
                .map(|d| d.0)
                .collect::<Vec<f64>>(),
            &self.average_cap,
        ))
    }

    /// Return temperature spread under test conditions, interpolated between design flow temps
    ///
    /// Arguments:
    /// * `flow_temp` - flow temp in K
    pub fn temp_spread_test_conditions(
        &self,
        design_flow_temp_op_cond: f64,
    ) -> Result<f64, BelowAbsoluteZeroError> {
        if self.design_flow_temps.len() == 1 {
            return Ok(self.temp_spread_test_conditions[0]);
        }

        let flow_temp = kelvin_to_celsius(design_flow_temp_op_cond)?;
        Ok(np_interp(
            flow_temp,
            &self
                .design_flow_temps
                .iter()
                .map(|d| d.0)
                .collect::<Vec<f64>>(),
            &self.temp_spread_test_conditions,
        ))
    }

    fn find_test_record_index(
        &self,
        test_condition: &str,
        design_flow_temp: f64,
    ) -> anyhow::Result<usize> {
        if test_condition == "cld" {
            return Ok(0);
        }

        self.test_data[&OrderedFloat(design_flow_temp)]
            .iter()
            .position(|test_record| test_record.test_letter.to_string() == test_condition)
            .ok_or_else(|| anyhow!("Test record index not found."))
    }

    /// Return value at specified test condition, interpolated between design flow temps
    fn data_at_test_condition(
        &self,
        data_item_name: DatumItem,
        test_condition: &str,
        design_flow_temp_op_cond: f64,
    ) -> Result<f64, BelowAbsoluteZeroError> {
        if self.design_flow_temps.len() == 1 {
            let idx = self
                .find_test_record_index(test_condition, self.design_flow_temps[0].0)
                .expect("Expected a test condition to be found in the test data for heat pumps");
            return Ok(self.test_data[&self.design_flow_temps[0]][idx].data_item(&data_item_name));
        }

        let data_list = &self
            .design_flow_temps
            .iter()
            .map(|design_flow_temp| {
                let idx = self
                    .find_test_record_index(test_condition, design_flow_temp.0)
                    .expect(
                        "Expected a test condition to be found in the test data for heat pumps",
                    );
                self.test_data[design_flow_temp][idx].data_item(&data_item_name)
            })
            .collect::<Vec<_>>();

        let flow_temp = kelvin_to_celsius(design_flow_temp_op_cond)?;
        Ok(np_interp(
            flow_temp,
            &self
                .design_flow_temps
                .iter()
                .map(|d| d.0)
                .collect::<Vec<f64>>(),
            data_list,
        ))
    }

    /// Return Carnot CoP at specified test condition (A, B, C, D, F or cld),
    /// interpolated between design flow temps
    fn carnot_cop_at_test_condition(
        &self,
        test_condition: &str,
        design_flow_temp_op_cond: f64,
    ) -> anyhow::Result<f64> {
        Ok(self.data_at_test_condition(
            DatumItem::CarnotCop,
            test_condition,
            design_flow_temp_op_cond,
        )?)
    }

    /// Return outlet temp, in Kelvin, at specified test condition (A, B, C, D,
    /// F or cld), interpolated between design flow temps.
    fn outlet_temp_at_test_condition(
        &self,
        test_condition: &str,
        design_flow_temp_op_cond: f64,
    ) -> Result<f64, BelowAbsoluteZeroError> {
        celsius_to_kelvin(self.data_at_test_condition(
            DatumItem::TempOutlet,
            test_condition,
            design_flow_temp_op_cond,
        )?)
    }

    /// Return source temp, in Kelvin, at specified test condition (A, B, C, D,
    /// F or cld), interpolated between design flow temps.
    fn source_temp_at_test_condition(
        &self,
        test_condition: &str,
        design_flow_temp_op_cond: f64,
    ) -> Result<f64, BelowAbsoluteZeroError> {
        celsius_to_kelvin(self.data_at_test_condition(
            DatumItem::TempSource,
            test_condition,
            design_flow_temp_op_cond,
        )?)
    }

    #[cfg(test)]
    /// Return capacity, in kW, at specified test condition (A, B, C, D, F or
    /// cld), interpolated between design flow temps.
    fn capacity_at_test_condition(
        &self,
        test_condition: &str,
        design_flow_temp_op_cond: f64,
    ) -> Result<f64, BelowAbsoluteZeroError> {
        self.data_at_test_condition(
            DatumItem::Capacity,
            test_condition,
            design_flow_temp_op_cond,
        )
    }

    /// Return load ratio at operating conditions
    fn load_ratio_at_operating_conditions(
        &self,
        temp_output: f64,
        temp_source: f64,
        carnot_cop_op_cond: f64,
        design_flow_temp: f64,
    ) -> anyhow::Result<f64> {
        let design_flow_temp = celsius_to_kelvin(design_flow_temp)?;
        let temp_output_cld = self.outlet_temp_at_test_condition("cld", design_flow_temp)?;
        let temp_source_cld = self.source_temp_at_test_condition("cld", design_flow_temp)?;
        let carnot_cop_cld = self.carnot_cop_at_test_condition("cld", design_flow_temp)?;

        Ok((carnot_cop_op_cond / carnot_cop_cld)
            * (temp_output_cld * temp_source / (temp_output * temp_source_cld)).powf(N_EXER))
    }

    /// Return test results either side of operating conditions.
    ///
    /// This function returns 4 results:
    /// - Exergy load ratio below operating conditions
    /// - Exergy load ratio above operating conditions
    /// - Exergy efficiency below operating conditions
    /// - Exergy efficiency above operating conditions
    ///
    /// Arguments:
    /// * `design_flow_temp` - flow temperature, in Kelvin
    /// * `exergy_lr_op_cond` - exergy load ratio at operating conditions
    fn lr_eff_either_side_of_op_cond(
        &self,
        design_flow_temp: &OrderedFloat<f64>,
        exergy_lr_op_cond: f64,
    ) -> anyhow::Result<(f64, f64, f64, f64)> {
        // For each design flow temperature, find load ratios in test data
        // either side of load ratio calculated for operating conditions.
        let design_flow_temp_data = &self.test_data[design_flow_temp];
        // Find the first load ratio in the test data that is greater than
        // or equal to than the load ratio at operating conditions - this
        // and the previous load ratio are the values either side of
        // operating conditions.
        let idx = match design_flow_temp_data
            .iter()
            .position(|test_record| test_record.theoretical_load_ratio > exergy_lr_op_cond)
        {
            Some(index) => {
                if index == 0 {
                    // Use second lowest (list index 1) and lowest
                    1
                } else {
                    index
                }
            }
            None => {
                // Use the highest (list index -1) and second highest
                design_flow_temp_data.len() - 1
            }
        };

        // Look up correct load ratio and efficiency based on the idx found above
        let lr_below = design_flow_temp_data[idx - 1].theoretical_load_ratio;
        let lr_above = design_flow_temp_data[idx].theoretical_load_ratio;
        let eff_below = design_flow_temp_data[idx - 1].exergetic_eff;
        let eff_above = design_flow_temp_data[idx].exergetic_eff;

        Ok((lr_below, lr_above, eff_below, eff_above))
    }

    fn exer_eff_op_cond(
        &self,
        temp_output: f64,
        temp_source: f64,
        carnot_cop_op_cond: f64,
        design_flow_temp_op_cond: f64,
    ) -> anyhow::Result<f64> {
        let mut exer_eff_op_cond_list: Vec<f64> = Default::default();

        for design_flow_temp in &self.design_flow_temps {
            // Get exergy load ratio at operating conditions
            // Get exergy load ratio and exergy efficiency at test conditions
            // above and below operating conditions
            let lr_op_cond = self.load_ratio_at_operating_conditions(
                temp_output,
                temp_source,
                carnot_cop_op_cond,
                **design_flow_temp,
            )?;
            let (lr_below, lr_above, eff_below, eff_above) =
                self.lr_eff_either_side_of_op_cond(design_flow_temp, lr_op_cond)?;

            // CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.4
            // Get exergy efficiency by interpolating between figures above and
            // below operating conditions
            let exer_eff_op_cond = if lr_below == lr_above {
                eff_below
            } else {
                eff_below
                    + (eff_below - eff_above) * (lr_op_cond - lr_below) / (lr_below - lr_above)
            };

            exer_eff_op_cond_list.push(exer_eff_op_cond);
        }

        let flow_temp = kelvin_to_celsius(design_flow_temp_op_cond)?;
        let design_temps_for_interp = &self
            .design_flow_temps
            .iter()
            .map(|temp| temp.0)
            .collect::<Vec<_>>();

        Ok(np_interp(
            flow_temp,
            design_temps_for_interp,
            &exer_eff_op_cond_list,
        ))
    }

    /// Calculate CoP at operating conditions when heat pump is not air-source
    ///
    /// Arguments:
    /// * `temp_diff_limit_low` - minimum temperature difference between source and sink
    /// * `temp_ext` - external temperature, in Celsius. Need to use Celsius for temp_ext_C
    ///                because regression coeffs were calculated using temperature in Celsius.
    /// * `temp_source` - source temperature, in Kelvin
    /// * `temp_output` - output temperature, in Kelvin
    /// * `design_flow_temp_op_cond` - design flow temperature of the heat pump as installed
    fn cop_op_cond_if_not_air_source(
        &self,
        temp_diff_limit_low: f64,
        temp_ext_c: f64,
        temp_source: f64,
        temp_output: f64,
        design_flow_temp_op_cond: f64,
    ) -> anyhow::Result<f64> {
        // For each design flow temperature, calculate CoP at operating conditions
        // Note: Loop over sorted list of design flow temps and then index into
        //        self.test_data, rather than looping over self.test_data,
        //        which is unsorted and therefore may populate the lists in the
        //        wrong order.
        let cop_op_cond = self
            .design_flow_temps
            .iter()
            .map(|design_flow_temp| {
                let design_flow_temp_data = &self.test_data[design_flow_temp];
                // Get the source and outlet temperatures from the coldest test record
                let temp_outlet_cld = celsius_to_kelvin(design_flow_temp_data[0].temp_outlet)?;
                let temp_source_cld = celsius_to_kelvin(design_flow_temp_data[0].temp_source)?;

                Ok((self.regression_coeffs[design_flow_temp][0]
                    + self.regression_coeffs[design_flow_temp][1] * temp_ext_c
                    + self.regression_coeffs[design_flow_temp][2] * temp_ext_c.powi(2))
                    * temp_output
                    * (temp_outlet_cld - temp_source_cld)
                    / (temp_outlet_cld * max_of_2(temp_output - temp_source, temp_diff_limit_low)))
            })
            .collect::<anyhow::Result<Vec<_>>>()?;

        if self.design_flow_temps.len() == 1 {
            return Ok(cop_op_cond[0]);
        }

        // Interpolate between the values found for the different design flow temperatures
        let flow_temp = kelvin_to_celsius(design_flow_temp_op_cond)?;
        Ok(np_interp(
            flow_temp,
            &self
                .design_flow_temps
                .iter()
                .map(|temp| temp.0)
                .collect::<Vec<_>>(),
            &cop_op_cond,
        ))
    }

    /// Calculate thermal capacity at operating conditions when heat pump is not air-source
    ///
    /// Arguments:
    /// * `temp_source` - source temperature, in Kelvin
    /// * `temp_output` - output temperature, in Kelvin
    /// * `mod_ctrl` - boolean specifying whether or not the heat has controls
    ///                    capable of varying the output (as opposed to just on/off
    ///                    control)
    /// * `design_flow_temp_op_cond` - design flow temperature of the heat pump as installed
    fn capacity_op_cond_var_flow_or_source_temp(
        &self,
        temp_output: f64,
        temp_source: f64,
        mod_ctrl: bool,
        design_flow_temp_op_cond: f64,
    ) -> Result<f64, BelowAbsoluteZeroError> {
        // In eqns below, method uses condition A rather than coldest. From
        // CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.4:
        // The Temperature Operation Limit (TOL) is defined in EN14825 as
        // "the lowest outdoor temperature at which the unit can still
        // deliver heating capacity and is declared by the manufacturer.
        // Below this temperature the heat pump will not be able to
        // deliver any heating capacity."
        // The weather data used within this calculation method does not
        // feature a source temperature at or below the "TOL" test
        // temperature (which is -7C to -10C). Therefore, test data at
        // the TOL test condition is not used (Test condition "A" at -7C
        // is sufficient).
        // TODO (from Python) The above implies that the TOL test temperature data may
        //       be needed if we change the weather data from that used in
        //       DAHPSE for SAP 2012/10.2
        let therm_cap_op_cond =
            self.design_flow_temps.iter().map(|design_flow_temp| {
                let design_flow_temp_data = &self.test_data[design_flow_temp];
                // Get the source and outlet temperatures from the coldest test record
                let temp_outlet_cld = celsius_to_kelvin(design_flow_temp_data[0].temp_outlet)?;
                let temp_source_cld = celsius_to_kelvin(design_flow_temp_data[0].temp_source)?;
                // Get the thermal capacity from the coldest test record
                let thermal_capacity_cld = design_flow_temp_data[0].capacity;

                Ok(if mod_ctrl {
                    thermal_capacity_cld * ((temp_outlet_cld * temp_source) / (temp_output * temp_source_cld)).powf(N_EXER)
                } else {
                    let d_idx = self.find_test_record_index("D", design_flow_temp.0).expect("Expected to find a test record with the condition 'D' within heat pump test data");
                    // Get the source and outlet temperatures for test condition D
                    let temp_outlet_d = celsius_to_kelvin(design_flow_temp_data[d_idx].temp_outlet)?;
                    let temp_source_d = celsius_to_kelvin(design_flow_temp_data[d_idx].temp_source)?;
                    // Get the thermal capacity for test condition D
                    let thermal_capacity_d = design_flow_temp_data[d_idx].capacity;

                    let temp_diff_cld = temp_outlet_cld - temp_source_cld;
                    let temp_diff_d = temp_outlet_d - temp_source_d;
                    let temp_diff_op_cond = temp_output - temp_source;

                    thermal_capacity_cld + (thermal_capacity_d - thermal_capacity_cld) * ((temp_diff_cld - temp_diff_op_cond) / (temp_diff_cld - temp_diff_d))
                })
            }).collect::<Result<Vec<f64>, BelowAbsoluteZeroError>>()?;

        // Interpolate between the values found for the different design flow temperatures
        let flow_temp = kelvin_to_celsius(design_flow_temp_op_cond)?;
        Ok(np_interp(
            flow_temp,
            &self
                .design_flow_temps
                .iter()
                .map(|temp| temp.0)
                .collect::<Vec<_>>(),
            &therm_cap_op_cond,
        ))
    }

    /// Calculate temperature spread correction factor
    ///
    /// Arguments:
    /// * `temp_source` - source temperature, in Kelvin
    /// * `temp_output` - output temperature, in Kelvin
    /// * `temp_diff_evaporator`
    ///     - average temperature difference between heat transfer medium and
    ///               refrigerant in evaporator, in deg C or Kelvin
    /// * `temp_diff_condenser`
    ///     - average temperature difference between heat transfer medium and
    ///               refrigerant in condenser, in deg C or Kelvin
    /// * `temp_spread_emitter`
    ///     - temperature spread on condenser side in operation due to design
    ///               of heat emission system
    /// * `design_flow_temp_op_cond`
    ///     - design flow temperature of the heat pump as installed
    fn temp_spread_correction(
        &self,
        temp_source: f64,
        temp_output: f64,
        temp_diff_evaporator: f64,
        temp_diff_condenser: f64,
        temp_spread_emitter: f64,
        design_flow_temp_op_cond: f64,
    ) -> Result<f64, BelowAbsoluteZeroError> {
        let temp_spread_correction_list = self
            .design_flow_temps
            .iter()
            .enumerate()
            .map(|(i, _design_flow_temp)| {
                let temp_spread_test_cond = self.temp_spread_test_conditions[i];
                1. - ((temp_spread_test_cond - temp_spread_emitter) / 2.)
                    / (temp_output - temp_spread_test_cond / 2. + temp_diff_condenser - temp_source
                        + temp_diff_evaporator)
            })
            .collect::<Vec<_>>();

        let flow_temp = kelvin_to_celsius(design_flow_temp_op_cond)?;
        Ok(np_interp(
            flow_temp,
            &self
                .design_flow_temps
                .iter()
                .map(|d| d.0)
                .collect::<Vec<f64>>(),
            &temp_spread_correction_list,
        ))
    }
}

/// The list average_cap will be in the same order as the
/// corresponding elements in self.__design_flow_temps. This behaviour
/// is relied upon elsewhere.
fn average_capacity(
    design_flow_temps: &[OrderedFloat<f64>],
    test_data: &HashMap<OrderedFloat<f64>, Vec<HeatPumpTestDatum>>,
) -> Vec<f64> {
    design_flow_temps
        .iter()
        .map(|design_flow_temp| {
            test_data[design_flow_temp]
                .iter()
                .filter(|datum| datum.test_letter.is_non_bivalent())
                .map(|datum| datum.capacity)
                .sum::<f64>()
                / TEST_LETTERS_NON_BIVALENT.len() as f64
        })
        .collect()
}

const TEMP_SPREAD_AT_TEST_CONDITIONS: [(f64, f64); 5] =
    [(20., 5.0), (35., 5.0), (45., 6.0), (55., 8.0), (65., 10.0)];

/// List temp spread at test conditions for the design flow temps in the test data
fn init_temp_spread_test_conditions(
    design_flow_temps: &[OrderedFloat<f64>],
) -> anyhow::Result<Vec<f64>> {
    design_flow_temps
        .iter()
        .map(|design_flow_temp| {
            TEMP_SPREAD_AT_TEST_CONDITIONS
                .iter()
                .find_map(|(temp_flow, temp_spread)| {
                    if design_flow_temp == temp_flow {
                        Some(*temp_spread)
                    } else {
                        None
                    }
                })
                .ok_or(anyhow!("Unexpected design flow temperature encountered"))
        })
        .collect()
}

fn init_regression_coeffs(
    design_flow_temps: &Vec<OrderedFloat<f64>>,
    test_data: &HashMap<OrderedFloat<f64>, Vec<HeatPumpTestDatum>>,
) -> anyhow::Result<HashMap<OrderedFloat<f64>, Vec<f64>>> {
    let mut regression_coeffs: HashMap<OrderedFloat<f64>, Vec<f64>> = Default::default();
    for design_flow_temp in design_flow_temps {
        let temp_test_list: Vec<f64> = test_data[design_flow_temp]
            .iter()
            .map(|datum| datum.temp_test)
            .collect();
        let cop_list: Vec<f64> = test_data[design_flow_temp]
            .iter()
            .map(|datum| datum.cop)
            .collect();
        regression_coeffs.insert(
            *design_flow_temp,
            polyfit(&temp_test_list, &cop_list, 2).map_err(|err| anyhow!(err))?,
        );
    }

    Ok(regression_coeffs)
}

impl PartialEq for HeatPumpTestDatum {
    fn eq(&self, other: &Self) -> bool {
        self.temp_test == other.temp_test && self.design_flow_temp == other.design_flow_temp
    }
}

// Values for time constant from BS EN 15316-4-2:2017 Table 13
const TIME_CONSTANT_WATER: f64 = 1560.;

/// An object to represent a water heating service provided by a heat pump to e.g. a cylinder.
///
/// This object contains the parts of the heat pump calculation that are
/// specific to providing hot water.
#[derive(Clone, Debug)]
pub struct HeatPumpServiceWater {
    heat_pump: Arc<Mutex<HeatPump>>,
    service_name: String,
    control: Arc<Control>,
    control_min: Arc<Control>,
    control_max: Arc<Control>,
    temp_limit_upper_in_k: f64,
    cold_feed: Arc<WaterSourceWithTemperature>,
    hybrid_boiler_service: Option<Arc<Mutex<BoilerServiceWaterRegular>>>,
}

impl HeatPumpServiceWater {
    pub(crate) fn new(
        heat_pump: Arc<Mutex<HeatPump>>,
        service_name: String,
        temp_limit_upper_in_c: f64,
        cold_feed: Arc<WaterSourceWithTemperature>,
        control_min: Arc<Control>, // TODO in Python 1.0.0a1 this is TimeControl
        control_max: Arc<Control>, // TODO in Python 1.0.0a1 this is TimeControl
        boiler_service_water_regular: Option<Arc<Mutex<BoilerServiceWaterRegular>>>,
    ) -> Self {
        let control = control_min.clone();
        Self {
            heat_pump,
            service_name,
            control,
            control_min,
            control_max,
            temp_limit_upper_in_k: celsius_to_kelvin(temp_limit_upper_in_c).expect(
                "Upper temp limit for heat pump is never expected to be below absolute zero.",
            ),
            cold_feed,
            hybrid_boiler_service: boiler_service_water_regular,
        }
    }

    pub fn is_on(&self, simtime: SimulationTimeIteration) -> bool {
        self.control.is_on(simtime)
    }

    /// Return water heating setpoint (not necessarily temperature)
    pub(crate) fn setpnt(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, Option<f64>)> {
        // TODO review - this check may be able to be removed in future if we validate control earlier
        match (&self.control_min.as_ref(), &self.control_max.as_ref()) {
            (
                Control::CombinationTime { .. } | Control::SetpointTime { .. },
                Control::CombinationTime { .. } | Control::SetpointTime { .. },
            ) => Ok((
                self.control_min.setpnt(&simulation_time_iteration),
                self.control_max.setpnt(&simulation_time_iteration),
            )),
            _ => bail!(
                "Expected control_min and control_max to be combination or setpoint time controls"
            ),
        }
    }

    const SERVICE_TYPE: HeatingServiceType = HeatingServiceType::DomesticHotWaterRegular;

    /// Calculate the maximum energy output of the HP, accounting for time
    /// spent on higher-priority services
    pub fn energy_output_max(
        &self,
        temp_flow: f64,
        temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Option<BufferTankEmittersDataWithResult>)> {
        // In the Python, temp_return has been updated to optional in 0.34 but update may be erroneous so leaving as non-optional for now
        // In 1.0.0a1 it looks like making temp_return optional is more intentional/correct due to a check of whether it's None at the beginning of this method before it's used:
        //         if not self.is_on() or temp_return is None:
        //             return 0.0
        // but wherever it's called, a float(non-optional) is passed in so leaving non-optional TODO review in later migrations
        let temp_return_k = celsius_to_kelvin(temp_return)?;
        if !self.is_on(simulation_time_iteration) {
            return Ok((0.0, None));
        }
        let temp_flow_k = celsius_to_kelvin(temp_flow)?;
        let design_flow_temp_op_cond = temp_flow_k;

        self.heat_pump.lock().energy_output_max(
            temp_flow_k,
            temp_return_k,
            design_flow_temp_op_cond,
            self.hybrid_boiler_service
                .as_ref()
                .map(|service| HybridBoilerService::Regular(service.clone())),
            Some(Self::SERVICE_TYPE),
            None,
            None,
            None,
            None,
            simulation_time_iteration,
        )
    }

    /// Demand energy (in kWh) from the heat pump
    pub fn demand_energy(
        &self,
        energy_demand: f64,
        temp_flow: Option<f64>,
        temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // In the Python, temp_return has been updated to optional in 0.34 but update may be erroneous so leaving as non-optional for now
        // (it's passed from here through methods to run_demand_energy_calc, where it is treated as non-optional and will cause an error if none
        // above is still the case in 1.0.0a1 TODO review in later migrations
        let service_on = self.is_on(simulation_time_iteration);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        let temp_return_k = celsius_to_kelvin(temp_return)?;

        if temp_flow.is_none() && energy_demand != 0. {
            bail!("temp_flow is None and energy_demand is not 0")
        };

        let temp_flow_k = if let Some(temp_flow) = temp_flow {
            Some(celsius_to_kelvin(temp_flow)?)
        } else {
            None
        };
        let design_flow_temp_op_cond = temp_flow_k.ok_or_else(|| anyhow!("temp_flow_k is None"))?;

        // TODO (from Python) Arbitrary volume used here for reference cold water temperature
        let list_temp_vol = self
            .cold_feed
            .get_temp_cold_water(1., simulation_time_iteration)?;
        let sum_t_by_v: f64 = list_temp_vol.iter().map(|(t, v)| t * v).sum();
        let sum_v: f64 = list_temp_vol.iter().map(|(_, v)| v).sum();
        let temp_cold_water = celsius_to_kelvin(sum_t_by_v / sum_v)?;

        self.heat_pump.lock().demand_energy(
            &self.service_name,
            &Self::SERVICE_TYPE,
            energy_demand,
            temp_flow_k,
            temp_return_k,
            self.temp_limit_upper_in_k,
            design_flow_temp_op_cond,
            TIME_CONSTANT_WATER,
            service_on,
            simulation_time_iteration,
            None,
            Some(temp_cold_water),
            None,
            self.hybrid_boiler_service
                .as_ref()
                .map(|service| HybridBoilerService::Regular(service.clone())),
            None,
            None,
        )
    }
}

// Values for time constant from BS EN 15316-4-2:2017 Table 13
const TIME_CONSTANT_SPACE_UFH: f64 = 1370.;
const TIME_CONSTANT_SPACE_FAN_COILS: f64 = 360.;
const TIME_CONSTANT_SPACE_WARM_AIR: f64 = 120.;

#[derive(Clone, Debug)]
pub struct HeatPumpServiceSpace {
    heat_pump: Arc<Mutex<HeatPump>>,
    service_name: String,
    emitter_type: HeatPumpEmitterType,
    control: Arc<Control>,
    temp_limit_upper_in_k: f64,
    temp_diff_emit_design: f64,
    design_flow_temp_op_cond: f64,
    hybrid_boiler_service: Option<Arc<Mutex<BoilerServiceSpace>>>,
    volume_heated: f64,
}

/// An object to represent a space heating service provided by a heat pump to e.g. radiators.
///
/// This object contains the parts of the heat pump calculation that are
/// specific to providing space heating.
impl HeatPumpServiceSpace {
    const SERVICE_TYPE: HeatingServiceType = HeatingServiceType::Space;

    /// Arguments:
    /// * `heat_pump` - reference to the HeatPump object providing the service
    /// * `service_name` - name of the service demanding energy from the heat pump
    /// * `emitter_type` - type of emitters that this space heating service is serving
    /// * `temp_limit_upper` - upper operating limit for temperature, in deg C
    /// * `temp_diff_emit_design` - design temperature difference across the emitters, in deg C or K
    /// * `design_flow_temp_op_cond` - design flow temperature of the heat pump as installed
    /// * `control` - reference to a control object which must implement is_on() and setpnt() funcs
    /// * `volume_heated` - volume of zones heated (required for exhaust air HPs only), in m3
    /// * `boiler_service_space` - reference to the BoilerServiceWaterSpace object that is backing up or supplementing the heat pump service
    pub(crate) fn new(
        heat_pump: Arc<Mutex<HeatPump>>,
        service_name: String,
        emitter_type: HeatPumpEmitterType,
        temp_limit_upper_in_c: f64,
        temp_diff_emit_design: f64,
        design_flow_temp_op_cond: f64,
        control: Arc<Control>, // TODO in Python 1.0.0a1 this is TimeControl | None (not making this optional in Rust 1.0.0a1 as there will always be a control)
        volume_heated: f64,
        boiler_service_space: Option<Arc<Mutex<BoilerServiceSpace>>>,
    ) -> Self {
        Self {
            heat_pump,
            service_name,
            emitter_type,
            control,
            temp_limit_upper_in_k: celsius_to_kelvin(temp_limit_upper_in_c)
                .expect("Upper temp limit for heat pump never expected to be below absolute zero"),
            temp_diff_emit_design,
            design_flow_temp_op_cond: celsius_to_kelvin(design_flow_temp_op_cond)
                .expect("design_flow_temp_op_cond expected to be above absolute zero"),
            hybrid_boiler_service: boiler_service_space,
            volume_heated,
        }
    }

    pub fn is_on(&self, simtime: SimulationTimeIteration) -> bool {
        self.control.is_on(simtime)
    }

    pub fn temp_setpnt(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<Option<f64>> {
        // TODO review - this check may be able to be removed in future if we validate control earlier
        match self.control.as_ref() {
            Control::CombinationTime { .. } | Control::SetpointTime { .. } => {
                Ok(self.control.setpnt(simulation_time_iteration))
            }
            _ => bail!("Expected control to be combination or setpoint time control"),
        }
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<Option<bool>> {
        // TODO review - this check may be able to be removed in future if we validate control earlier
        match self.control.as_ref() {
            Control::CombinationTime { .. } | Control::SetpointTime { .. } => {
                Ok(self.control.in_required_period(simulation_time_iteration))
            }
            _ => bail!("Expected control to be combination or setpoint time control"),
        }
    }

    /// Calculate the maximum energy output of the HP, accounting for time
    /// spent on higher-priority services
    pub fn energy_output_max(
        &self,
        temp_output: f64,
        temp_return_feed: f64,
        time_start: Option<f64>,
        emitters_data_for_buffer_tank: Option<BufferTankEmittersData>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Option<BufferTankEmittersDataWithResult>)> {
        let time_start = time_start.unwrap_or(0.);
        if !self.is_on(simulation_time_iteration) {
            return Ok((0.0, None));
        }

        let temp_output_k = celsius_to_kelvin(temp_output)?;
        let mut heat_pump = self.heat_pump.lock();
        let source_type = heat_pump.source_type;
        heat_pump.energy_output_max(
            temp_output_k,
            temp_return_feed,
            self.design_flow_temp_op_cond,
            self.hybrid_boiler_service
                .as_ref()
                .map(|service| HybridBoilerService::Space(service.clone())),
            Some(Self::SERVICE_TYPE),
            Some(TempSpreadCorrectionArg::Callable(
                self.temp_spread_correction_fn(source_type),
            )),
            Some(time_start),
            emitters_data_for_buffer_tank,
            Some(self.service_name.as_str()),
            simulation_time_iteration,
        )
    }

    /// Demand energy (in kWh) from the heat pump
    ///
    /// Arguments:
    /// * `energy_demand` - space heating energy demand, in kWh
    /// * `temp_flow` - flow temperature for emitters, in deg C
    /// * `temp_return` - return temperature for emitters, in deg C
    /// * `update_heat_source_state` - if false the heat pump state does not change
    pub fn demand_energy(
        &self,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        time_start: Option<f64>,
        emitters_data_for_buffer_tank: Option<BufferTankEmittersDataWithResult>,
        update_heat_source_state: Option<bool>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);
        let service_on = self.is_on(simulation_time_iteration);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        let mut heat_pump = self.heat_pump.lock();
        let time_constant_for_service = match self.emitter_type {
            HeatPumpEmitterType::RadiatorsUfh => TIME_CONSTANT_SPACE_UFH,
            HeatPumpEmitterType::FanCoils => TIME_CONSTANT_SPACE_FAN_COILS,
            HeatPumpEmitterType::WarmAir => TIME_CONSTANT_SPACE_WARM_AIR,
        };
        let source_type = heat_pump.source_type;
        heat_pump.demand_energy(
            &self.service_name,
            &Self::SERVICE_TYPE,
            energy_demand,
            Some(celsius_to_kelvin(temp_flow)?),
            celsius_to_kelvin(temp_return)?,
            self.temp_limit_upper_in_k,
            self.design_flow_temp_op_cond,
            time_constant_for_service,
            service_on,
            simulation_time_iteration,
            Some(TempSpreadCorrectionArg::Callable(
                self.temp_spread_correction_fn(source_type),
            )),
            None,
            Some(time_start),
            self.hybrid_boiler_service
                .as_ref()
                .map(|service| HybridBoilerService::Space(service.clone())),
            emitters_data_for_buffer_tank,
            Some(update_heat_source_state),
        )
    }

    /// Return the cumulative running time and throughput factor (exhaust air HPs only)
    pub fn running_time_throughput_factor(
        &self,
        space_heat_running_time_cumulative: f64,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        time_start: Option<f64>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        let time_start = time_start.unwrap_or(0.);
        let service_on = self.is_on(simulation_time_iteration);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };
        let mut heat_pump = self.heat_pump.lock();
        let time_constant_for_service = match self.emitter_type {
            HeatPumpEmitterType::RadiatorsUfh => TIME_CONSTANT_SPACE_UFH,
            HeatPumpEmitterType::FanCoils => TIME_CONSTANT_SPACE_FAN_COILS,
            HeatPumpEmitterType::WarmAir => TIME_CONSTANT_SPACE_WARM_AIR,
        };
        let source_type = heat_pump.source_type;

        heat_pump.running_time_throughput_factor(
            space_heat_running_time_cumulative,
            &self.service_name,
            &HeatingServiceType::Space,
            energy_demand,
            celsius_to_kelvin(temp_flow)?,
            celsius_to_kelvin(temp_return)?,
            self.temp_limit_upper_in_k,
            self.design_flow_temp_op_cond,
            time_constant_for_service,
            service_on,
            self.volume_heated,
            simulation_time_iteration,
            Some(TempSpreadCorrectionArg::Callable(
                self.temp_spread_correction_fn(source_type),
            )),
            Some(time_start),
        )
    }

    fn temp_spread_correction_fn(
        &self,
        source_type: HeatPumpSourceType,
    ) -> Box<dyn Fn(f64, f64, f64, HeatPumpTestData) -> f64> {
        // Average temperature difference between heat transfer medium and
        // refrigerant in condenser
        let temp_diff_condenser = 5.0;

        // Average temperature difference between heat transfer medium and
        // refrigerant in evaporator
        let temp_diff_evaporator = match source_type {
            t if t.source_fluid_is_air() => 15.0,
            t if t.source_fluid_is_water() => 10.0,
            _ => panic!("impossible heat pump source type encountered"),
        };

        let temp_diff_emit_design = self.temp_diff_emit_design;

        Box::new(
            move |temp_output,
                  temp_source,
                  design_flow_temp_op_cond,
                  test_data: HeatPumpTestData| {
                test_data.temp_spread_correction(
                    temp_source,
                    temp_output,
                    temp_diff_evaporator,
                    temp_diff_condenser,
                    temp_diff_emit_design,
                    design_flow_temp_op_cond,
                ).expect("Did not expect temp spread correction function to encounter an illegal temperature")
            },
        )
    }
}

/// An object to represent a warm air space heating service provided by a heat pump.
///
///    This object contains the parts of the heat pump calculation that are
///    specific to providing space heating via warm air.
#[derive(Clone, Debug)]
pub struct HeatPumpServiceSpaceWarmAir {
    heat_pump: Arc<Mutex<HeatPump>>,
    service_name: String,
    control: Arc<Control>,
    emitter_type: HeatPumpEmitterType,
    temp_limit_upper_in_k: f64,
    temp_diff_emit_design: f64,
    design_flow_temp_op_cond: f64,
    frac_convective: f64,
    temp_flow: f64,
    temp_return: f64,
    volume_heated: f64,
}

impl HeatPumpServiceSpaceWarmAir {
    pub(crate) fn new(
        heat_pump: Arc<Mutex<HeatPump>>,
        service_name: &str,
        temp_diff_emit_design: f64,
        design_flow_temp_op_cond: f64,
        control: Arc<Control>, // TODO in Python 1.0.0a1 this is TimeControl | None (not making this optional in Rust 1.0.0a1 as there will always be a control)
        temp_flow: f64,
        frac_convective: f64,
        volume_heated: f64,
    ) -> Self {
        let temp_limit_upper_in_c = temp_flow;

        Self {
            heat_pump,
            service_name: service_name.into(),
            emitter_type: HeatPumpEmitterType::WarmAir,
            control,
            temp_limit_upper_in_k: celsius_to_kelvin(temp_limit_upper_in_c).expect(
                "Upper limit given for heat pump never expected to be below absolute zero.",
            ),
            temp_diff_emit_design,
            design_flow_temp_op_cond,
            frac_convective,
            temp_flow,
            temp_return: temp_flow,
            volume_heated,
        }
    }

    pub fn is_on(&self, simtime: SimulationTimeIteration) -> bool {
        self.control.is_on(simtime)
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<Option<bool>> {
        // TODO review - this check may be able to be removed in future if we validate control earlier
        match self.control.as_ref() {
            Control::CombinationTime { .. } | Control::SetpointTime { .. } => {
                Ok(self.control.in_required_period(simulation_time_iteration))
            }
            _ => bail!("Expected control to be combination or setpoint time control"),
        }
    }

    pub fn temp_setpnt(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<Option<f64>> {
        // TODO review - this check may be able to be removed in future if we validate control earlier
        match self.control.as_ref() {
            Control::CombinationTime { .. } | Control::SetpointTime { .. } => {
                Ok(self.control.setpnt(simulation_time_iteration))
            }
            _ => bail!("Expected control to be combination or setpoint time control"),
        }
    }

    pub(crate) fn energy_output_min(&self) -> f64 {
        0.
    }

    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let temp_flow = self.temp_flow;
        let temp_return = self.temp_return;

        let service_on = self.is_on(simulation_time_iteration);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        let mut pump = self.heat_pump.lock();
        let time_constant_for_service = match self.emitter_type {
            HeatPumpEmitterType::RadiatorsUfh => TIME_CONSTANT_SPACE_UFH,
            HeatPumpEmitterType::FanCoils => TIME_CONSTANT_SPACE_FAN_COILS,
            HeatPumpEmitterType::WarmAir => TIME_CONSTANT_SPACE_WARM_AIR,
        };

        let source_type = pump.source_type;

        pump.demand_energy(
            &self.service_name,
            &HeatingServiceType::Space,
            energy_demand,
            Some(celsius_to_kelvin(temp_flow)?),
            celsius_to_kelvin(temp_return)?,
            self.temp_limit_upper_in_k,
            self.design_flow_temp_op_cond,
            time_constant_for_service,
            service_on,
            simulation_time_iteration,
            Some(TempSpreadCorrectionArg::Callable(
                self.temp_spread_correction_fn(source_type),
            )),
            None,
            None,
            None,
            None,
            None,
        )
    }

    // TODO: will this ever be called if only relevant for exhaust air heat pump and
    // given that create_service_space_heating_warm_air bails if heat pump sink_type is not air?
    /// Return the cumulative running time and throughput factor (exhaust air HPs only)
    /// Arguments:
    /// energy_demand -- in kWh
    /// space_heat_running_time_cumulative
    /// -- running time spent on higher-priority space heating services
    pub fn running_time_throughput_factor(
        &self,
        energy_demand: f64,
        space_heat_running_time_cumulative: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        let temp_flow = self.temp_flow;
        let temp_return = self.temp_return;

        let service_on = self.is_on(simulation_time_iteration);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        let time_constant_for_service = match self.emitter_type {
            HeatPumpEmitterType::RadiatorsUfh => TIME_CONSTANT_SPACE_UFH,
            HeatPumpEmitterType::FanCoils => TIME_CONSTANT_SPACE_FAN_COILS,
            HeatPumpEmitterType::WarmAir => TIME_CONSTANT_SPACE_WARM_AIR,
        };
        let mut pump = self.heat_pump.lock();
        let source_type = pump.source_type;

        pump.running_time_throughput_factor(
            space_heat_running_time_cumulative,
            &self.service_name,
            &HeatingServiceType::Space,
            energy_demand,
            celsius_to_kelvin(temp_flow)?,
            celsius_to_kelvin(temp_return)?,
            self.temp_limit_upper_in_k,
            self.design_flow_temp_op_cond,
            time_constant_for_service,
            service_on,
            self.volume_heated,
            simulation_time_iteration,
            Some(TempSpreadCorrectionArg::Callable(
                self.temp_spread_correction_fn(source_type),
            )),
            None,
        )
    }

    /// yes this is copy-pasted from the SpaceService
    fn temp_spread_correction_fn(
        &self,
        source_type: HeatPumpSourceType,
    ) -> Box<dyn Fn(f64, f64, f64, HeatPumpTestData) -> f64> {
        // Average temperature difference between heat transfer medium and
        // refrigerant in condenser
        let temp_diff_condenser = 5.0;

        // Average temperature difference between heat transfer medium and
        // refrigerant in evaporator
        let temp_diff_evaporator = match source_type {
            t if t.source_fluid_is_air() => 15.0,
            t if t.source_fluid_is_water() => 10.0,
            _ => panic!("impossible heat pump source type encountered"),
        };

        let temp_diff_emit_design = self.temp_diff_emit_design;

        Box::new(
            move |temp_output,
                  temp_source,
                  design_flow_temp_op_cond,
                  test_data: HeatPumpTestData| {
                test_data.temp_spread_correction(
                    temp_source,
                    temp_output,
                    temp_diff_evaporator,
                    temp_diff_condenser,
                    temp_diff_emit_design,
                    design_flow_temp_op_cond,
                ).expect("Temp spread correction never expects to encounter temperature below absolute zero")
            },
        )
    }

    pub fn frac_convective(&self) -> f64 {
        self.frac_convective
    }
}

/// From CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.3:
/// A minimum temperature difference of 6K between the source and sink
/// temperature is applied to prevent very high Carnot COPs entering the
/// calculation. This only arises when the temperature difference and heating
/// load is small and is unlikely to affect the calculated SPF.
const HEAT_PUMP_TEMP_DIFF_LIMIT_LOW: f64 = 6.0; // Kelvin

/// Fraction of the energy input dedicated to auxiliaries when on
// TODO (from Python) This is always zero for electric heat pumps, but if we want to deal
//                    with non-electric heat pumps then this will need to be altered.
const HEAT_PUMP_F_AUX: f64 = 0.0;

/// A type to represent an electric heat pump
#[derive(Clone, Derivative)]
#[derivative(Debug)]
pub struct HeatPump {
    // energy supply
    pub(crate) energy_supply: Arc<RwLock<EnergySupply>>,
    energy_supply_connections: HashMap<String, Arc<EnergySupplyConnection>>,
    energy_supply_connection_aux: Arc<EnergySupplyConnection>,
    simulation_timestep: f64,
    external_conditions: Arc<ExternalConditions>,
    // litres/ second
    service_results: Arc<RwLock<Vec<ServiceResult>>>,
    total_time_running_current_timestep_full_load: f64,
    time_running_continuous: f64,
    source_type: HeatPumpSourceType,
    sink_type: HeatPumpSinkType,
    backup_ctrl: HeatPumpBackupControlType,
    time_delay_backup: Option<f64>,
    modulating_ctrl: bool,
    time_constant_onoff_operation: f64,
    temp_return_feed_max: Option<f64>,
    temp_lower_op_limit: f64,
    temp_diff_flow_return_min: f64,
    var_flow_temp_ctrl_during_test: bool,
    power_heating_circ_pump: f64,
    power_heating_warm_air_fan: f64,
    power_source_circ_pump: f64,
    power_standby: f64,
    power_crankcase_heater_mode: f64,
    power_off_mode: f64,
    power_max_backup: f64,
    temp_distribution_heat_network: Option<f64>,
    // energy supply for heat networks
    energy_supply_heat_source: Option<Arc<RwLock<EnergySupply>>>,
    boiler: Option<Arc<RwLock<Boiler>>>,
    cost_schedule_hybrid_hp: Option<BoilerCostScheduleHybrid>,
    cost_schedule_metadata: Option<CostScheduleMetadata>,
    volume_heated_all_services: Option<f64>,
    eahp_mixed_max_temp: Option<f64>,
    eahp_mixed_min_temp: Option<f64>,
    ext_air_ratio: Option<f64>,
    buffer_tank: Option<BufferTank>,
    #[derivative(Debug = "ignore")]
    temp_internal_air_fn: TempInternalAirFn,
    energy_supply_heat_source_connections: HashMap<String, Arc<EnergySupplyConnection>>,
    overvent_ratio: f64,
    test_data: HeatPumpTestData,
    temp_min_modulation_rate_low: Option<f64>,
    min_modulation_rate_low: Option<f64>,
    temp_min_modulation_rate_high: Option<f64>,
    min_modulation_rate_55: Option<f64>,
    detailed_results: Option<Vec<Arc<RwLock<Vec<ServiceResult>>>>>,
    #[cfg(test)]
    canned_whether_delay: Option<bool>,
}

impl HeatPump {
    /// Construct a HeatPump object
    ///
    /// Arguments:
    /// * `heat_pump_input` - An input object containing heat pump details
    ///     NB. must be heat pump variant or this will panic
    /// * `energy_supply` - reference to EnergySupply object
    /// * `energy_supply_conn_name_auxiliary`
    ///             - name to be used for EnergySupplyConnection object for auxiliary energy
    /// * `external_conditions` - reference to ExternalConditions object
    /// * `number_of_zones` - number of zones in the project
    /// * `throughput_exhaust_air` - throughput (litres / second) of exhaust air
    /// * `energy_supply_heat_source` - reference to EnergySupply object representing heat network
    ///                        (for HPs that use heat network as heat source) / other heat source
    /// * `output_detailed_results` - if true, save detailed results from each timestep
    ///                                   for later reporting
    /// * `boiler` - reference to Boiler object representing boiler (for hybrid heat pumps)
    /// * `cost_schedule_hybrid_hp` - cost schedules for heat pump and boiler
    ///
    /// Other variables:
    ///        energy_supply_connections
    ///            -- dictionary with service name strings as keys and corresponding
    ///               EnergySupplyConnection objects as values
    ///        energy_supply_connection_aux -- EnergySupplyConnection object for auxiliary energy
    ///        test_data -- HeatPumpTestData object
    pub(crate) fn new(
        heat_pump_input: &HeatSourceWetDetails,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_conn_name_auxiliary: &str,
        simulation_timestep: f64,
        external_conditions: Arc<ExternalConditions>,
        number_of_zones: usize,
        throughput_exhaust_air: Option<f64>,
        energy_supply_heat_source: Option<Arc<RwLock<EnergySupply>>>,
        output_detailed_results: bool,
        boiler: Option<Arc<RwLock<Boiler>>>,
        cost_schedule_hybrid_hp: Option<BoilerCostScheduleHybrid>,
        temp_internal_air_fn: TempInternalAirFn,
    ) -> anyhow::Result<Self> {
        let energy_supply_connections = Default::default();
        let energy_supply_connection_aux = Arc::new(EnergySupply::connection(
            energy_supply.clone(),
            energy_supply_conn_name_auxiliary,
        )?);

        let service_results = Default::default();
        let total_time_running_current_timestep_full_load = Default::default();
        let time_running_continuous = Default::default();

        let (
            source_type,
            sink_type,
            backup_ctrl,
            time_delay_backup,
            modulating_ctrl,
            time_constant_onoff_operation,
            temp_return_feed_max,
            temp_lower_op_limit,
            temp_diff_flow_return_min,
            var_flow_temp_ctrl_during_test,
            power_heating_circ_pump,
            power_heating_warm_air_fan,
            power_source_circ_pump,
            power_standby,
            power_crankcase_heater_mode,
            power_off_mode,
            power_max_backup,
            temp_distribution_heat_network,
            test_data,
            min_modulation_rate_20,
            min_modulation_rate_35,
            min_modulation_rate_55,
            eahp_mixed_max_temp,
            eahp_mixed_min_temp,
            buffer_tank,
        ) = match heat_pump_input {
            HeatSourceWetDetails::HeatPump {
                source_type,
                sink_type,
                backup_control_type,
                time_delay_backup,
                modulating_control,
                time_constant_onoff_operation,
                temp_return_feed_max,
                temp_lower_operating_limit,
                min_temp_diff_flow_return_for_hp_to_operate,
                var_flow_temp_ctrl_during_test,
                power_heating_circ_pump,
                power_heating_warm_air_fan,
                power_source_circ_pump,
                power_standby,
                power_crankcase_heater,
                power_off,
                power_max_backup,
                temp_distribution_heat_network,
                test_data_en14825: test_data,
                min_modulation_rate_20,
                min_modulation_rate_35,
                min_modulation_rate_55,
                eahp_mixed_max_temp,
                eahp_mixed_min_temp,
                buffer_tank,
                ..
            } => (
                *source_type,
                *sink_type,
                *backup_control_type,
                *time_delay_backup,
                *modulating_control,
                *time_constant_onoff_operation,
                *temp_return_feed_max,
                celsius_to_kelvin(*temp_lower_operating_limit)?,
                *min_temp_diff_flow_return_for_hp_to_operate,
                *var_flow_temp_ctrl_during_test,
                *power_heating_circ_pump,
                *(power_heating_warm_air_fan.as_ref().unwrap_or(&0.)),
                *power_source_circ_pump,
                *power_standby,
                *power_crankcase_heater,
                *power_off,
                *power_max_backup,
                temp_distribution_heat_network,
                test_data.clone(),
                min_modulation_rate_20,
                min_modulation_rate_35,
                min_modulation_rate_55,
                *eahp_mixed_max_temp,
                *eahp_mixed_min_temp,
                buffer_tank,
            ),
            _ => {
                bail!("Expected heat source details passed to heat pump to be of heat pump variant")
            }
        };
        let (time_delay_backup, power_max_backup) = match backup_ctrl {
            HeatPumpBackupControlType::None => (None, Some(0.)),
            _ => (
                time_delay_backup,
                match boiler {
                    None => power_max_backup,
                    Some(_) => Some(0.),
                },
            ),
        };

        let cost_schedule_metadata = cost_schedule_hybrid_hp
            .as_ref()
            .map(|cost_schedule| anyhow::Ok(cost_schedule.try_into()?))
            .transpose()?;

        let temp_return_feed_max = match sink_type {
            HeatPumpSinkType::Air => None,
            _ => temp_return_feed_max.map(celsius_to_kelvin).transpose()?,
        };
        let temp_distribution_heat_network = match source_type {
            HeatPumpSourceType::HeatNetwork => *temp_distribution_heat_network,
            _ => None,
        };

        // Exhaust air HP requires different/additional initialisation, which is implemented here
        let (overvent_ratio, volume_heated_all_services, test_data_after_interpolation) =
            match throughput_exhaust_air {
                Some(throughput_exhaust_air) if source_type.is_exhaust_air() => {
                    let (lowest_air_flow_rate_in_test_data, test_data_after_interpolation) =
                        interpolate_exhaust_air_heat_pump_test_data(
                            throughput_exhaust_air,
                            &test_data,
                            source_type,
                        )?;
                    (
                        max_of_2(
                            1.0,
                            lowest_air_flow_rate_in_test_data / throughput_exhaust_air,
                        ),
                        Some(f64::default()),
                        test_data_after_interpolation,
                    )
                }
                _ => (1.0, None, test_data.into_iter().map(Into::into).collect()),
            };

        // TODO (from Python) For now, disable exhaust air heat pump when conditions are out of
        //                    range of test data. Decision to be made on whether to allow this
        //                    and if so how to handle it in the calculation
        if overvent_ratio > 1.0 {
            bail!("Exhaust air heat pump: flow rate for associated ventilation system is lower than flow rate in any of the test data provided.");
        }

        // in the Python there is a check here about an air flow rate being present - decided
        // any validation should really happen upstream of here

        let (eahp_mixed_max_temp, eahp_mixed_min_temp, ext_air_ratio) =
            if source_type == HeatPumpSourceType::ExhaustAirMixed {
                let ext_air_ratio_list: Vec<_> = test_data_after_interpolation
                    .iter()
                    .filter_map(|datum| datum.eahp_mixed_ext_air_ratio)
                    .collect();
                let same_ext_air_ratio = ext_air_ratio_list
                    .iter()
                    .map(|f| OrderedFloat(*f))
                    .unique()
                    .count();
                if same_ext_air_ratio > 1 {
                    bail!("More than one unique external air ratio entered.");
                }
                let ext_air_ratio: f64 =
                    ext_air_ratio_list.iter().sum::<f64>() / ext_air_ratio_list.len() as f64;
                (
                    eahp_mixed_max_temp,
                    eahp_mixed_min_temp,
                    Some(ext_air_ratio),
                )
            } else {
                (None, None, None)
            };

        let test_data = HeatPumpTestData::new(test_data_after_interpolation)?;

        let (
            temp_min_modulation_rate_low,
            min_modulation_rate_low,
            temp_min_modulation_rate_high,
            min_modulation_rate_55,
        ) = if modulating_ctrl {
            let (temp_min_modulation_rate_low, min_modulation_rate_low) = match sink_type {
                HeatPumpSinkType::Air => (
                    Some(celsius_to_kelvin(20.)?),
                    Some(
                        min_modulation_rate_20
                            .ok_or(anyhow!("expected a min_modulation_rate_20 provided"))?,
                    ),
                ),
                HeatPumpSinkType::Water | HeatPumpSinkType::Glycol25 => (
                    Some(celsius_to_kelvin(35.)?),
                    Some(
                        min_modulation_rate_35
                            .ok_or(anyhow!("expected a min_modulation_rate_35 provided"))?,
                    ),
                ),
            };
            let (temp_min_modulation_rate_high, min_modulation_rate_55) =
                if test_data.design_flow_temps.contains(&OrderedFloat(55.)) {
                    (
                        Some(celsius_to_kelvin(55.)?),
                        Some(
                            min_modulation_rate_55
                                .ok_or(anyhow!("expected a min_modulation_rate_55 provided"))?,
                        ),
                    )
                } else {
                    (None, None)
                };

            (
                temp_min_modulation_rate_low,
                min_modulation_rate_low,
                temp_min_modulation_rate_high,
                min_modulation_rate_55,
            )
        } else {
            (None, None, None, None)
        };

        let detailed_results = if output_detailed_results {
            Some(vec![])
        } else {
            None
        };

        let buffer_tank = buffer_tank.as_ref().map(|tank| {
            BufferTank::new(
                tank.daily_losses,
                tank.volume,
                tank.pump_fixed_flow_rate,
                tank.pump_power_at_flow_rate,
                number_of_zones,
                simulation_timestep,
                temp_internal_air_fn(),
                *WATER,
                Some(output_detailed_results),
            )
        });

        Ok(Self {
            energy_supply,
            energy_supply_connections,
            energy_supply_connection_aux,
            simulation_timestep,
            external_conditions,
            service_results,
            total_time_running_current_timestep_full_load,
            time_running_continuous,
            source_type,
            sink_type,
            backup_ctrl,
            time_delay_backup,
            modulating_ctrl,
            time_constant_onoff_operation,
            temp_return_feed_max,
            temp_lower_op_limit,
            temp_diff_flow_return_min,
            var_flow_temp_ctrl_during_test,
            power_heating_circ_pump: power_heating_circ_pump.unwrap_or(0.),
            power_heating_warm_air_fan,
            power_source_circ_pump,
            power_standby,
            power_crankcase_heater_mode,
            power_off_mode,
            power_max_backup: power_max_backup
                .expect("A power_max_backup value was expected but not provided."),
            temp_distribution_heat_network,
            energy_supply_heat_source,
            boiler,
            cost_schedule_hybrid_hp,
            cost_schedule_metadata,
            volume_heated_all_services,
            eahp_mixed_max_temp,
            eahp_mixed_min_temp,
            ext_air_ratio,
            buffer_tank,
            temp_internal_air_fn,
            energy_supply_heat_source_connections: Default::default(),
            overvent_ratio,
            test_data,
            temp_min_modulation_rate_low,
            min_modulation_rate_low,
            temp_min_modulation_rate_high,
            min_modulation_rate_55,
            detailed_results,
            #[cfg(test)]
            canned_whether_delay: None,
        })
    }

    #[cfg(test)]
    pub(crate) fn backup_heater_max_power(&self) -> f64 {
        self.power_max_backup
    }

    pub fn source_is_exhaust_air(&self) -> bool {
        self.source_type.is_exhaust_air()
    }

    pub fn buffer_int_gains(&self) -> f64 {
        if let Some(buffer_tank) = self.buffer_tank.as_ref() {
            buffer_tank.internal_gains()
        } else {
            0.0
        }
    }

    fn create_service_connection(
        heat_pump: Arc<Mutex<Self>>,
        service_name: &str,
    ) -> Result<(), anyhow::Error> {
        // Check that the service_name is not already registered
        if heat_pump
            .lock()
            .energy_supply_connections
            .contains_key(service_name)
        {
            bail!("Error: service name already used: {service_name}");
        }

        // Set up energy supply connection for this service
        {
            let mut heat_pump = heat_pump.lock();
            let connection = Arc::new(EnergySupply::connection(
                heat_pump.energy_supply.clone(),
                service_name,
            )?);
            heat_pump
                .energy_supply_connections
                .insert(service_name.into(), connection);
        }

        let energy_supply_heat_source = heat_pump.lock().energy_supply_heat_source.clone();

        if let Some(energy_supply_heat_source) = energy_supply_heat_source {
            heat_pump
                .lock()
                .energy_supply_heat_source_connections
                .insert(
                    service_name.into(),
                    Arc::new(EnergySupply::connection(
                        energy_supply_heat_source.clone(),
                        service_name,
                    )?),
                );
        }

        Ok(())
    }

    pub(crate) fn create_service_hot_water_combi(
        &self,
        boiler_data: HotWaterSourceDetails,
        service_name: &str,
        temp_hot_water: f64,
        cold_feed: WaterSourceWithTemperature,
    ) -> anyhow::Result<BoilerServiceWaterCombi> {
        if let Some(boiler) = self.boiler.as_ref() {
            Boiler::create_service_hot_water_combi(
                boiler.clone(),
                boiler_data,
                service_name,
                temp_hot_water,
                cold_feed,
            )
            .map_err(|err| anyhow!(err))
        } else {
            bail!("Error: Missing Boiler object for hybrid heat pump. Non-hybrid heat pumps cannot be used as an instantaneous hot water system.")
        }
    }

    pub(crate) fn create_service_hot_water(
        heat_pump: Arc<Mutex<Self>>,
        service_name: &str,
        temp_limit_upper_in_c: f64,
        cold_feed: Arc<WaterSourceWithTemperature>,
        control_min: Arc<Control>, // TODO in Python 1.0.0a1 this is SetpointTimeControl | CombinationTimeControl
        control_max: Arc<Control>, // TODO in Python 1.0.0a1 this is SetpointTimeControl | CombinationTimeControl
    ) -> anyhow::Result<HeatPumpServiceWater> {
        Self::create_service_connection(heat_pump.clone(), service_name)?;
        let boiler_service = heat_pump
            .lock()
            .boiler
            .as_ref()
            .map(|boiler| {
                anyhow::Ok(Arc::new(Mutex::new(
                    Boiler::create_service_hot_water_regular(
                        boiler.clone(),
                        service_name,
                        control_min.clone(),
                        control_max.clone(),
                    )?,
                )))
            })
            .transpose()?;

        Ok(HeatPumpServiceWater::new(
            heat_pump,
            service_name.into(),
            temp_limit_upper_in_c,
            cold_feed,
            control_min,
            control_max,
            boiler_service,
        ))
    }

    pub(crate) fn create_service_space_heating(
        heat_pump: Arc<Mutex<Self>>,
        service_name: &str,
        emitter_type: HeatPumpEmitterType,
        temp_limit_upper_in_c: f64,
        temp_diff_emit_design: f64,
        design_flow_temp_op_cond: f64,
        control: Arc<Control>, // TODO in Python 1.0.0a1 this is SetpointTimeControl | CombinationTimeControl | None (not making this optional in Rust 1.0.0a1 as there will always be a control)
        volume_heated: f64,
    ) -> HeatPumpServiceSpace {
        let boiler_service = heat_pump.lock().boiler.as_ref().map(|boiler| {
            Arc::new(Mutex::new(Boiler::create_service_space_heating(
                boiler.clone(),
                service_name,
                control.clone(),
            )))
        });
        if heat_pump.lock().source_is_exhaust_air() {
            if let Some(volume_heated_all_services) =
                heat_pump.lock().volume_heated_all_services.as_mut()
            {
                *volume_heated_all_services += volume_heated;
            }
        }
        Self::create_service_connection(heat_pump.clone(), service_name).unwrap();
        HeatPumpServiceSpace::new(
            heat_pump,
            service_name.into(),
            emitter_type,
            temp_limit_upper_in_c,
            temp_diff_emit_design,
            design_flow_temp_op_cond,
            control,
            volume_heated,
            boiler_service,
        )
    }

    pub(crate) fn create_service_space_heating_warm_air(
        heat_pump: Arc<Mutex<Self>>,
        service_name: &str,
        control: Arc<Control>, // TODO in Python 1.0.0a1 this is SetpointTimeControl | CombinationTimeControl | None (not making this optional in Rust 1.0.0a1 as there will always be a control)
        frac_convective: f64,
        volume_heated: f64,
    ) -> anyhow::Result<HeatPumpServiceSpaceWarmAir> {
        {
            let mut heat_pump = heat_pump.lock();
            if heat_pump.sink_type != HeatPumpSinkType::Air {
                bail!("Warm air space heating service requires heat pump with sink type Air");
            }
            // following is erroneous logic that is following upstream for now
            if heat_pump.boiler.is_some() {
                // TODO (from Python) More evidence is required before warm air space heating systems can work with hybrid heat pumps
                bail!("Cannot handle hybrid warm air heat pumps - calculation not implemented");
            }
            if heat_pump.source_is_exhaust_air() {
                if let Some(volume_heated_all_services) =
                    heat_pump.volume_heated_all_services.as_mut()
                {
                    *volume_heated_all_services += volume_heated;
                }
            }
        }

        // Use low temperature test data for space heating - set flow temp such
        // that it matches the one used in the test
        let temp_flow = heat_pump.lock().test_data.design_flow_temps[0].0;
        let design_flow_temp_op_cond = temp_flow;
        // Design temperature difference across the emitters, in deg C or K
        let temp_diff_emit_design = max_of_2(
            temp_flow / 7.,
            heat_pump
                .lock()
                .test_data
                .temp_spread_test_conditions(design_flow_temp_op_cond)?,
        );

        Self::create_service_connection(heat_pump.clone(), service_name)?;
        Ok(HeatPumpServiceSpaceWarmAir::new(
            heat_pump,
            service_name,
            temp_diff_emit_design,
            design_flow_temp_op_cond,
            control,
            temp_flow,
            frac_convective,
            volume_heated,
        ))
    }

    /// Get source temp according to rules in CALCM-01 - DAHPSE - V2.0_DRAFT13, 3.1.1
    fn get_temp_source(&self, simulation_time_iteration: SimulationTimeIteration) -> f64 {
        let temp_source = match self.source_type {
            HeatPumpSourceType::Ground => {
                let temp_ext = self.external_conditions.air_temp(&simulation_time_iteration);
                max_of_2(
                    0.,
                    min_of_2(8., temp_ext * 0.25806 + 2.8387),
                )
            }
            HeatPumpSourceType::OutsideAir => self.external_conditions.air_temp(&simulation_time_iteration),
            HeatPumpSourceType::ExhaustAirMEV => (self.temp_internal_air_fn)(),
            HeatPumpSourceType::ExhaustAirMVHR => (self.temp_internal_air_fn)(),
            HeatPumpSourceType::ExhaustAirMixed => {
                let temp_ext = self.external_conditions.air_temp(&simulation_time_iteration);
                let temp_int = (self.temp_internal_air_fn)();
                let ext_air_ratio = self.ext_air_ratio.expect("An ext_air_ratio is expected when the heat pump source type is ExhaustAirMixed");
                let temp_mixed = ext_air_ratio * temp_ext + (1. - ext_air_ratio) * temp_int;
                if temp_ext > self.eahp_mixed_max_temp.expect("An eahp_mixed_max_temp is expected when the heat pump source type is ExhaustAirMixed") || temp_mixed < self.eahp_mixed_min_temp.expect("An eahp_mixed_min_temp is expected when the heat pump source type is ExhaustAirMixed") {
                    temp_int
                } else {
                    temp_mixed
                }
            }
            HeatPumpSourceType::WaterGround => {
                // Where water is extracted from the ground and re-injected into the
                // ground or discharged at the surface, the surface temperature is assumed
                // to be constant and equal to the annual average air temperature.
                self.external_conditions.air_temp_annual().expect("An annual air temp is expected to be accessible.")
            }
            HeatPumpSourceType::WaterSurface => {
                // Where water is extracted from surface water, such as rivers and lakes,
                // it is assumed that this extraction does not substantively affect the
                // average temperature of the water volume, thus it must have a sufficient
                // thermal capacity. The source temperature is taken as the monthly air
                // temperature.
                self.external_conditions.air_temp_monthly(simulation_time_iteration.current_month_start_end_hours())
            }
            HeatPumpSourceType::HeatNetwork => self.temp_distribution_heat_network.expect("expected a temperature to be available for a heat network if this is the source type"),
        };

        celsius_to_kelvin(temp_source)
            .expect("In context, temperature is never expected to be less than absolute zero.")
    }

    /// Calculate the thermal capacity of the heat pump at operating conditions
    ///
    /// Based on CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.4
    fn thermal_capacity_op_cond(
        &self,
        temp_output: f64,
        temp_source: f64,
        design_flow_temp_op_cond: f64,
    ) -> Result<f64, BelowAbsoluteZeroError> {
        if !matches!(self.source_type, HeatPumpSourceType::OutsideAir)
            && !self.var_flow_temp_ctrl_during_test
        {
            self.test_data.average_capacity(design_flow_temp_op_cond)
        } else {
            self.test_data.capacity_op_cond_var_flow_or_source_temp(
                temp_output,
                temp_source,
                self.modulating_ctrl,
                design_flow_temp_op_cond,
            )
        }
    }

    fn backup_energy_output_max(
        &self,
        temp_output: f64,
        temp_return_feed: f64,
        time_available: f64,
        time_start: f64,
        hybrid_boiler_service: Option<HybridBoilerService>,
        simtime: SimulationTimeIteration,
    ) -> Result<f64, BelowAbsoluteZeroError> {
        let time_elapsed_hp = match self.backup_ctrl {
            HeatPumpBackupControlType::TopUp => None,
            HeatPumpBackupControlType::None | HeatPumpBackupControlType::Substitute => {
                Some(self.total_time_running_current_timestep_full_load)
            }
        };

        Ok(match hybrid_boiler_service {
            Some(HybridBoilerService::Regular(service_water_regular)) => {
                service_water_regular.lock().energy_output_max(
                    kelvin_to_celsius(temp_output)?,
                    kelvin_to_celsius(temp_return_feed)?,
                    time_elapsed_hp,
                    simtime,
                )
            }
            Some(HybridBoilerService::Space(service_space)) => {
                service_space.lock().energy_output_max(
                    kelvin_to_celsius(temp_output)?,
                    kelvin_to_celsius(temp_return_feed)?,
                    Some(time_start),
                    time_elapsed_hp,
                    simtime,
                )
            }
            None => self.power_max_backup * time_available,
        })
    }

    /// Calculate the maximum energy output of the HP, accounting for time
    /// spent on higher-priority services
    ///
    /// Note: Call via a HeatPumpService object, not directly.
    fn energy_output_max(
        &mut self,
        temp_output: f64,
        temp_return_feed: f64,
        design_flow_temp_op_cond: f64,
        hybrid_boiler_service: Option<HybridBoilerService>,
        service_type: Option<HeatingServiceType>,
        temp_spread_correction: Option<TempSpreadCorrectionArg>,
        time_start: Option<f64>,
        emitters_data_for_buffer_tank: Option<BufferTankEmittersData>,
        service_name: Option<&str>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Option<BufferTankEmittersDataWithResult>)> {
        let timestep = self.simulation_timestep;
        let time_start = time_start.unwrap_or(0.);
        let time_available = self.time_available(time_start, timestep, None);
        let temp_source = self.get_temp_source(simtime);
        let mut temp_output = temp_output;
        let temp_spread_correction =
            temp_spread_correction.unwrap_or(TempSpreadCorrectionArg::Float(1.0));

        // If there is buffer tank, energy output max will be affected by the temperature lift required
        // by the tank
        let mut heat_loss_buffer_kwh = 0.0;
        let emitters_data_for_buffer_tank = if let (Some(buffer_tank), Some(emitters_data)) =
            (self.buffer_tank.as_mut(), emitters_data_for_buffer_tank)
        {
            let buffer_tank_buffer_loss = buffer_tank.buffer_loss();
            let buffer_tank_results = buffer_tank.calc_buffer_tank(
                service_name.expect("A service name was expected to be available."),
                emitters_data,
            )?;
            let last_result = buffer_tank_results
                .iter()
                .last()
                .expect("At least one buffer tank result was expected.");

            // upstream Python checks for existence of flow_temp_increase_due_to_buffer and
            // heat_loss_buffer_kwh fields, but in Rust definition they are non-optional
            temp_output += last_result.flow_temp_increase_due_to_buffer;
            heat_loss_buffer_kwh = buffer_tank_buffer_loss + last_result.heat_loss_buffer_kwh;

            emitters_data.merge_result(last_result)
        } else {
            Default::default()
        };

        let hp_cost_effective = match (&self.cost_schedule_hybrid_hp, &self.boiler) {
            (Some(_), Some(boiler)) => {
                let cop_op_cond = self.cop_op_cond(
                    service_type,
                    temp_output,
                    temp_source,
                    temp_spread_correction,
                    design_flow_temp_op_cond,
                    simtime,
                )?;
                let energy_output_max_boiler = self.backup_energy_output_max(
                    temp_output,
                    temp_return_feed,
                    time_available,
                    time_start,
                    hybrid_boiler_service.clone(),
                    simtime,
                )?;
                let boiler_eff = boiler.read().calc_boiler_eff(
                    false,
                    kelvin_to_celsius(temp_return_feed)?,
                    energy_output_max_boiler,
                    Some(time_start),
                    Some(timestep),
                    simtime,
                )?;
                self.is_heat_pump_cost_effective(cop_op_cond, boiler_eff, simtime)
            }
            _ => true,
        };

        let energy_max = if !hp_cost_effective {
            self.backup_energy_output_max(
                temp_output,
                temp_return_feed,
                time_available,
                time_start,
                hybrid_boiler_service,
                simtime,
            )?
        } else {
            let outside_operating_limits = self.outside_operating_limits(temp_return_feed, simtime);
            let power_max_hp = if outside_operating_limits {
                0.0
            } else {
                self.thermal_capacity_op_cond(temp_output, temp_source, design_flow_temp_op_cond)?
            };

            match (
                self.backup_ctrl,
                self.backup_ctrl != HeatPumpBackupControlType::None
                    && self.backup_heater_delay_time_elapsed()
                    && outside_operating_limits,
            ) {
                (HeatPumpBackupControlType::None, _) | (_, false) => power_max_hp * time_available,
                (HeatPumpBackupControlType::TopUp, _) => {
                    let energy_max_backup = self.backup_energy_output_max(
                        temp_output,
                        temp_return_feed,
                        time_available,
                        time_start,
                        hybrid_boiler_service,
                        simtime,
                    )?;
                    (power_max_hp * time_available) + energy_max_backup
                }
                (HeatPumpBackupControlType::Substitute, _) => {
                    let energy_max_backup = self.backup_energy_output_max(
                        temp_output,
                        temp_return_feed,
                        time_available,
                        time_start,
                        hybrid_boiler_service,
                        simtime,
                    )?;
                    max_of_2(power_max_hp * time_available, energy_max_backup)
                }
            }
        };

        Ok((
            energy_max - heat_loss_buffer_kwh,
            if let (Some(_), emitters_data) =
                (self.buffer_tank.as_ref(), emitters_data_for_buffer_tank)
            {
                Some(emitters_data)
            } else {
                None
            },
        ))
    }

    /// Calculate CoP at operating conditions
    fn cop_op_cond(
        &self,
        _service_type: Option<HeatingServiceType>,
        temp_output: f64,
        temp_source: f64,
        temp_spread_correction: TempSpreadCorrectionArg,
        design_flow_temp_op_cond: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let temp_spread_correction_factor = match temp_spread_correction {
            TempSpreadCorrectionArg::Float(correction) => correction,
            TempSpreadCorrectionArg::Callable(callable) => callable(
                temp_output,
                temp_source,
                design_flow_temp_op_cond,
                self.test_data.clone(),
            ),
        };

        Ok(
            if !matches!(self.source_type, HeatPumpSourceType::OutsideAir)
                && !self.var_flow_temp_ctrl_during_test
            {
                temp_spread_correction_factor
                    * self.test_data.cop_op_cond_if_not_air_source(
                        HEAT_PUMP_TEMP_DIFF_LIMIT_LOW,
                        self.external_conditions
                            .air_temp(&simulation_time_iteration),
                        temp_source,
                        temp_output,
                        design_flow_temp_op_cond,
                    )?
            } else {
                let carnot_cop_op_cond = carnot_cop(
                    temp_source,
                    temp_output,
                    Some(HEAT_PUMP_TEMP_DIFF_LIMIT_LOW),
                );
                let exer_eff_op_cond = self.test_data.exer_eff_op_cond(
                    temp_output,
                    temp_source,
                    carnot_cop_op_cond,
                    design_flow_temp_op_cond,
                )?;

                // CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.5
                // Note: DAHPSE method document section 4.5.5 doesn't have
                // temp_spread_correction_factor in formula below. However, section 4.5.7
                // states that the correction factor is to be applied to the CoP.
                max_of_2(
                    1.0,
                    exer_eff_op_cond * carnot_cop_op_cond * temp_spread_correction_factor,
                )
            },
        )
    }

    /// Calculate energy output limited by upper temperature
    fn energy_output_limited(
        &self,
        energy_output_required: f64,
        temp_output: Option<f64>,
        temp_used_for_scaling: f64,
        temp_limit_upper: f64,
    ) -> f64 {
        if let (Some(temp_output), true) = (
            temp_output,
            temp_output.is_some_and(|temp_output| temp_output > temp_limit_upper),
        ) {
            if temp_output == temp_used_for_scaling {
                energy_output_required
            } else if (temp_limit_upper - temp_used_for_scaling) >= self.temp_diff_flow_return_min {
                energy_output_required * (temp_limit_upper - temp_used_for_scaling)
                    / (temp_output - temp_used_for_scaling)
            } else {
                0.
            }
        } else {
            energy_output_required
        }
    }

    #[cfg(test)]
    fn set_canned_whether_delay_time_elapsed(&mut self, whether: Option<bool>) {
        self.canned_whether_delay = whether;
    }

    /// Check if backup heater is available or still in delay period
    fn backup_heater_delay_time_elapsed(&self) -> bool {
        #[cfg(test)]
        {
            if let Some(whether) = self.canned_whether_delay {
                return whether;
            }
        }

        self.time_running_continuous
            >= self
                .time_delay_backup
                .expect("time delay backup was expected to be set")
    }

    /// Check if heat pump is outside operating limits
    fn outside_operating_limits(
        &self,
        temp_return_feed: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> bool {
        let temp_source = self.get_temp_source(simulation_time_iteration);
        let below_min_ext_temp = temp_source <= self.temp_lower_op_limit;

        let above_temp_return_feed_max = match self.sink_type {
            HeatPumpSinkType::Water | HeatPumpSinkType::Glycol25 => {
                temp_return_feed
                    > self
                        .temp_return_feed_max
                        .expect("a temp return feed max was expected to have been set")
            }
            HeatPumpSinkType::Air => false,
        };

        below_min_ext_temp || above_temp_return_feed_max
    }

    fn inadequate_capacity(
        &self,
        energy_output_required: f64,
        thermal_capacity_op_cond: f64,
        temp_output: f64,
        time_available: f64,
        time_start: f64,
        temp_return_feed: f64,
        hybrid_boiler_service: Option<HybridBoilerService>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<bool> {
        let timestep = self.simulation_timestep;

        // For top-up backup heater, use backup if delay time has elapsed.
        // For substitute backup heater, use backup if delay time has elapsed and
        // backup heater can provide more energy than heat pump. This assumption
        // is required to make the maximum energy output of the system
        // predictable before the demand is known.
        let energy_max_backup = self.backup_energy_output_max(
            temp_output,
            temp_return_feed,
            time_available,
            time_start,
            hybrid_boiler_service,
            simtime,
        )?;

        Ok(
            if (matches!(self.backup_ctrl, HeatPumpBackupControlType::TopUp)
                && self.backup_heater_delay_time_elapsed())
                || (matches!(self.backup_ctrl, HeatPumpBackupControlType::Substitute)
                    && self.backup_heater_delay_time_elapsed()
                    && energy_max_backup > (thermal_capacity_op_cond * time_available))
            {
                energy_output_required > thermal_capacity_op_cond * timestep
            } else {
                false
            },
        )
    }

    fn is_heat_pump_cost_effective(
        &self,
        cop_op_cond: f64,
        boiler_eff: f64,
        simtime: SimulationTimeIteration,
    ) -> bool {
        let schedule_metadata = self
            .cost_schedule_metadata
            .as_ref()
            .expect("Cost schedule information expected to be available.");
        let schedule_index = simtime.time_series_idx(
            schedule_metadata.start_day,
            schedule_metadata.time_series_step,
        );
        let cost_hp = schedule_metadata.hp[schedule_index];
        let cost_boiler = schedule_metadata.boiler[schedule_index];

        (cost_hp / cop_op_cond) <= (cost_boiler / boiler_eff)
    }

    /// Evaluate boolean conditions that may trigger backup heater
    fn use_backup_heater_only(
        &self,
        cop_op_cond: f64,
        energy_output_required: f64,
        thermal_capacity_op_cond: f64,
        temp_output: f64,
        time_available: f64,
        time_start: f64,
        temp_return_feed: f64,
        hybrid_boiler_service: Option<HybridBoilerService>,
        boiler_eff: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<bool> {
        let outside_operating_limits = self.outside_operating_limits(temp_return_feed, simtime);
        let inadequate_capacity = self.inadequate_capacity(
            energy_output_required,
            thermal_capacity_op_cond,
            temp_output,
            time_available,
            time_start,
            temp_return_feed,
            hybrid_boiler_service,
            simtime,
        )?;

        let hp_not_cost_effective = if self.cost_schedule_hybrid_hp.is_some() {
            !self.is_heat_pump_cost_effective(
                cop_op_cond,
                boiler_eff.expect("Boiler efficiency expected."),
                simtime,
            )
        } else {
            false
        };

        Ok(self.backup_ctrl != HeatPumpBackupControlType::None
            && (outside_operating_limits
                || (inadequate_capacity
                    && self.backup_ctrl == HeatPumpBackupControlType::Substitute)
                || hp_not_cost_effective))
    }

    /// Calculate time available for the current service
    fn time_available(
        &self,
        time_start: f64,
        timestep: f64,
        additional_time_unavailable: Option<f64>,
    ) -> f64 {
        // Assumes that time spent on other services is evenly spread throughout
        // the timestep so the adjustment for start time below is a proportional
        // reduction of the overall time available, not simply a subtraction
        let additional_time_unavailable = additional_time_unavailable.unwrap_or(0.);

        (timestep
            - self.total_time_running_current_timestep_full_load
            - additional_time_unavailable)
            * (1. - time_start / timestep)
    }

    /// Calculate energy required by heat pump to satisfy demand for the service indicated.
    ///
    /// Note: Call via the __demand_energy func, not directly.
    ///     This function should not save any results to member variables of
    ///     this class, because it may need to be run more than once (e.g. for
    ///     exhaust air heat pumps). Results should be returned to the
    ///     demand_energy function which calls this one and will save results
    ///     when appropriate.
    /// Note: The optional variable additional_time_unavailable is used for
    ///     calculating running time without needing to update any state - the
    ///     variable contains the time already committed to other services
    ///     where the running time has not been added to
    ///     self.total_time_running_current_timestep_full_load
    /// Note: The optional variable time_start is used to account for situations
    ///     where the service cannot start at the beginning of the timestep
    ///     (e.g. due to emitters having to cool down before the service can
    ///     run). This is then used to decrease the time available proportionally,
    ///     which assumes that the timing of other services is randomly
    ///     distributed, neither coinciding perfectly with the time that the
    ///     current service is running or the time it is not running. This
    ///     variable is relative to the beginning of the timestep
    fn run_demand_energy_calc(
        &mut self,
        service_name: &str,
        service_type: &HeatingServiceType,
        energy_output_required: f64,
        temp_output: Option<f64>,
        temp_return_feed: f64,
        temp_limit_upper: f64,
        design_flow_temp_op_cond: f64,
        time_constant_for_service: f64,
        service_on: bool,
        simtime: SimulationTimeIteration,
        temp_spread_correction: Option<TempSpreadCorrectionArg>,
        temp_used_for_scaling: Option<f64>,
        hybrid_boiler_service: Option<HybridBoilerService>,
        boiler_eff: Option<f64>,
        additional_time_unavailable: Option<f64>,
        time_start: Option<f64>,
        emitters_data_for_buffer_tank: Option<BufferTankEmittersDataWithResult>,
    ) -> anyhow::Result<HeatPumpEnergyCalculation> {
        let time_start = time_start.unwrap_or(0.);
        let (flow_temp_increase_due_to_buffer, power_buffer_tank_pump, heat_loss_buffer_kwh) =
            if let (Some(buffer_tank), Some(emitters_data)) =
                (self.buffer_tank.as_mut(), emitters_data_for_buffer_tank)
            {
                let buffer_tank_results = &emitters_data.result;
                let BufferTankServiceResult {
                    flow_temp_increase_due_to_buffer,
                    pump_power_at_flow_rate: power_buffer_tank_pump,
                    ..
                } = buffer_tank_results;
                let mut heat_loss_buffer_kwh =
                    buffer_tank.buffer_loss() + buffer_tank_results.heat_loss_buffer_kwh;
                if !service_on {
                    buffer_tank.update_buffer_loss(heat_loss_buffer_kwh);
                    heat_loss_buffer_kwh = 0.0;
                }
                (
                    *flow_temp_increase_due_to_buffer,
                    *power_buffer_tank_pump,
                    heat_loss_buffer_kwh,
                )
            } else {
                (0.0, 0.0, 0.0)
            };

        // Adding buffer tank losses to the energy required to be delivered by the HP
        let mut energy_output_required = energy_output_required;
        energy_output_required += heat_loss_buffer_kwh;

        let mut temp_output = temp_output;
        if temp_output.is_some() {
            temp_output = Some(temp_output.unwrap() + flow_temp_increase_due_to_buffer);
        }

        let temp_spread_correction =
            temp_spread_correction.unwrap_or(TempSpreadCorrectionArg::Float(1.0));
        let temp_used_for_scaling = temp_used_for_scaling.unwrap_or(temp_return_feed);
        let additional_time_unavailable = additional_time_unavailable.unwrap_or(0.0);

        let timestep = self.simulation_timestep;

        let energy_output_limited = self.energy_output_limited(
            energy_output_required,
            temp_output,
            temp_used_for_scaling,
            temp_limit_upper,
        );

        let temp_source = self.get_temp_source(simtime);
        // From here onwards, output temp to be used is subject to the upper limit
        if temp_output.is_some() {
            temp_output = Some(min_of_2(temp_output.unwrap(), temp_limit_upper));
        };

        // Get thermal capacity and CoP at operating conditions
        let (thermal_capacity_op_cond, cop_op_cond) = if let Some(temp_output) = temp_output {
            let thermal_capacity_op_cond =
                self.thermal_capacity_op_cond(temp_output, temp_source, design_flow_temp_op_cond)?;
            let cop_op_cond = self.cop_op_cond(
                Some(*service_type),
                temp_output,
                temp_source,
                temp_spread_correction,
                design_flow_temp_op_cond,
                simtime,
            )?;
            (Some(thermal_capacity_op_cond), Some(cop_op_cond))
        } else {
            (None, None)
        };

        // Calculate time available to run service
        // Note that different values are calculated here for calculations
        // related to capacity and calculations related to load ratio and
        // cycling.
        // For calculations related to capacity, we account for the time
        // required for emitters to cool to the maximum flow temperature at the
        // start of the timestep to account for the effect of needing to allow
        // the emitters to cool down on the heat pump's usable capacity. If
        // there is a requirement to allow emitters to cool down, then the heat
        // pump capacity will be limited by this requirement.
        // However, in reality the maximum flow temperature may be adjusted more
        // frequently than just once per modelling timestep, so the real effect
        // is likely to look more like an overall reduction in demand (and
        // therefore load ratio) across the whole timestep, and so when
        // calculating the load ratio we do not want to adjust the time
        // available as if there were a cool-down period only at the start of
        // the timestep, but instead consider the load ratio over the whole
        // timestep.
        let time_available_for_capacity_calc =
            self.time_available(time_start, timestep, Some(additional_time_unavailable));
        let time_available_for_load_ratio_calc =
            self.time_available(0., timestep, Some(additional_time_unavailable));

        let time_required = if let Some(thermal_capacity_op_cond) = thermal_capacity_op_cond {
            energy_output_limited / thermal_capacity_op_cond
        } else {
            0.
        };

        // TODO (from Python) Consider moving some of these checks earlier or to HeatPumpService
        //      classes. May be able to skip a lot of the calculation.
        let mut boiler_eff = boiler_eff;
        if let (Some(_), Some(boiler)) = (&hybrid_boiler_service, &self.boiler) {
            boiler_eff = Some(boiler.read().calc_boiler_eff(
                false,
                kelvin_to_celsius(temp_return_feed)?,
                energy_output_required,
                Some(time_start),
                Some(timestep),
                simtime,
            )?);
        }

        let use_backup_heater_only = match (temp_output, cop_op_cond, thermal_capacity_op_cond) {
            (Some(temp_output), Some(cop_op_cond), Some(thermal_capacity_op_cond)) => self
                .use_backup_heater_only(
                    cop_op_cond,
                    energy_output_required,
                    thermal_capacity_op_cond,
                    temp_output,
                    time_available_for_capacity_calc,
                    time_start,
                    temp_return_feed,
                    hybrid_boiler_service.clone(),
                    boiler_eff,
                    simtime,
                )?,
            _ => false,
        };

        // Calculate energy delivered by HP
        let (energy_delivered_hp, mut time_running_current_service_full_load) = match (
            self.compressor_is_running(service_on, use_backup_heater_only),
            thermal_capacity_op_cond,
        ) {
            (true, Some(thermal_capacity_op_cond)) => {
                // Backup heater not providing entire energy requirement
                let time_running_current_service_full_load =
                    min_of_2(time_required, time_available_for_capacity_calc);
                let energy_delivered_hp =
                    thermal_capacity_op_cond * time_running_current_service_full_load;

                (energy_delivered_hp, time_running_current_service_full_load)
            }
            _ => (0., 0.),
        };

        // Calculate energy delivered by backup heater
        let mut energy_delivered_backup = match (
            self.backup_ctrl,
            self.backup_ctrl != HeatPumpBackupControlType::None
                && !self.backup_heater_delay_time_elapsed()
                && !use_backup_heater_only
                || !service_on,
        ) {
            (HeatPumpBackupControlType::None, _) | (_, true) => 0.0,
            (HeatPumpBackupControlType::TopUp | HeatPumpBackupControlType::Substitute, _) => {
                let energy_max_backup = self.backup_energy_output_max(
                    temp_output.ok_or_else(|| anyhow!("Expected temp_output to be set"))?,
                    temp_return_feed,
                    time_available_for_capacity_calc,
                    time_start,
                    hybrid_boiler_service.clone(),
                    simtime,
                )?;
                let energy_delivered_backup = max_of_2(
                    min_of_2(
                        energy_max_backup,
                        energy_output_required - energy_delivered_hp,
                    ),
                    0.0,
                );
                // If all energy provided by backup heater, calculate running time consumed
                if use_backup_heater_only && self.power_max_backup > 0. {
                    time_running_current_service_full_load = min_of_2(
                        energy_delivered_backup / self.power_max_backup,
                        time_available_for_capacity_calc,
                    );
                }

                energy_delivered_backup
            }
        };

        // Calculate energy input to backup heater
        // TODO (from Python) Account for backup heater efficiency, or call another heating
        //      system object. For now, assume 100% efficiency
        let mut energy_input_backup = energy_delivered_backup;

        let energy_output_required_boiler = if hybrid_boiler_service.is_some() {
            let energy_output_required_boiler = energy_delivered_backup;
            energy_delivered_backup = 0.0;
            energy_input_backup = 0.0;
            energy_output_required_boiler
        } else {
            0.0
        };

        // Energy used by pumps
        let energy_source_circ_pump = if use_backup_heater_only {
            0.
        } else {
            time_running_current_service_full_load * self.power_source_circ_pump
        };

        let (energy_heating_warm_air_fan, energy_heating_circ_pump) = if *service_type
            == HeatingServiceType::Space
            && self.sink_type == HeatPumpSinkType::Air
        {
            // If warm air distribution add electricity for warm air fan
            (
                time_running_current_service_full_load * self.power_heating_warm_air_fan,
                0.,
            )
        } else {
            // if wet distribution add electricity for wet distribution
            (
                0.,
                time_running_current_service_full_load
                    * (self.power_heating_circ_pump + power_buffer_tank_pump),
            )
        };

        // Calculate total energy delivered and input
        let mut energy_delivered_total = energy_delivered_hp + energy_delivered_backup;

        if let Some(buffer_tank) = self.buffer_tank.as_mut() {
            if energy_delivered_total >= heat_loss_buffer_kwh {
                energy_delivered_total -= heat_loss_buffer_kwh;
                buffer_tank.update_buffer_loss(0.0);
            } else {
                let hp_covered_buffer_losses = energy_delivered_total;
                energy_delivered_total = 0.0;
                buffer_tank.update_buffer_loss(heat_loss_buffer_kwh - hp_covered_buffer_losses);
            }
        }

        Ok(HeatPumpEnergyCalculation {
            service_name: result_str(service_name),
            service_type: *service_type,
            service_on,
            energy_output_required,
            temp_output,
            temp_source,
            cop_op_cond,
            thermal_capacity_op_cond,
            time_running_full_load: time_running_current_service_full_load,
            time_running_part_load: None,
            time_available_for_capacity_calc,
            time_available_for_load_ratio_calc,
            time_constant_for_service,
            compressor_power_min_load: Default::default(),
            load_ratio_continuous_min: Default::default(),
            load_ratio: Default::default(),
            use_backup_heater_only,
            hp_operating_in_onoff_mode: Default::default(),
            energy_input_hp: Default::default(),
            energy_delivered_hp,
            energy_input_backup,
            energy_delivered_backup,
            energy_input_total: Default::default(),
            energy_delivered_total,
            energy_heating_circ_pump,
            energy_source_circ_pump,
            energy_output_required_boiler,
            energy_heating_warm_air_fan,
            energy_output_delivered_boiler: Default::default(),
        })
    }

    fn load_ratio_and_mode(
        &self,
        time_running_current_service: f64,
        time_available_for_current_service: f64,
        temp_output: f64,
    ) -> anyhow::Result<(f64, f64, bool)> {
        // Calculate load ratio
        let load_ratio = if is_close!(
            time_available_for_current_service,
            0.,
            rel_tol = 1e-09,
            abs_tol = 1e-10
        ) {
            if !is_close!(
                time_running_current_service,
                0.,
                rel_tol = 1e-09,
                abs_tol = 1e-10
            ) {
                bail!("Calculated time running is not zero despite no time being available")
            }
            0.
        } else {
            time_running_current_service / time_available_for_current_service
        };

        let load_ratio_continuous_min = if self.modulating_ctrl {
            let min_modulation_rate_low = self
                .min_modulation_rate_low
                .expect("A min_modulation_rate_low was expected to have been set");
            if self
                .test_data
                .design_flow_temps
                .contains(&OrderedFloat(55.))
            {
                let temp_min_modulation_rate_high = self
                    .temp_min_modulation_rate_high
                    .expect("temp_min_modulation_rate_high expected to have been provided");
                let temp_min_modulation_rate_low = self
                    .temp_min_modulation_rate_low
                    .expect("temp_min_modulation_rate_low expected to have been provided");
                np_interp(
                    temp_output,
                    &[temp_min_modulation_rate_low, temp_min_modulation_rate_high],
                    &[
                        min_modulation_rate_low,
                        self.min_modulation_rate_55
                            .expect("A min modulation rate was expected for a 55 air flow temp"),
                    ],
                )
            } else {
                min_modulation_rate_low
            }
        } else {
            // On/off heat pump cannot modulate below maximum power
            1.0
        };
        // Determine whether HP is operating in on/off mode
        let hp_operating_in_onoff_mode = load_ratio > 0.
            && !is_close!(load_ratio, 0., rel_tol = 1e-09, abs_tol = 1e-10)
            && load_ratio < load_ratio_continuous_min;

        Ok((
            load_ratio,
            load_ratio_continuous_min,
            hp_operating_in_onoff_mode,
        ))
    }

    fn compressor_is_running(&self, service_on: bool, use_backup_heater_only: bool) -> bool {
        service_on && !use_backup_heater_only
    }

    fn energy_input_compressor(
        &self,
        service_on: bool,
        use_backup_heater_only: bool,
        hp_operating_in_onoff_mode: bool,
        energy_delivered_hp: f64,
        energy_delivered_hp_aggregated: f64,
        thermal_capacity_op_cond: Option<f64>,
        cop_op_cond: f64,
        time_available_for_current_service: f64,
        load_ratio: f64,
        load_ratio_continuous_min: f64,
        time_constant_for_service: f64,
        _service_type: HeatingServiceType,
    ) -> (f64, f64) {
        let Some(thermal_capacity_op_cond) = thermal_capacity_op_cond else {
            return (0.0, 0.0);
        };

        let compressor_power_full_load = thermal_capacity_op_cond / cop_op_cond;

        // CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.10, step 1:
        let compressor_power_min_load = compressor_power_full_load * load_ratio_continuous_min;

        let energy_input_hp = if !self.compressor_is_running(service_on, use_backup_heater_only) {
            0.
        } else if hp_operating_in_onoff_mode {
            // BS EN 15316-4-2:2017 equation 28
            let power_used_due_to_inertia_effects = compressor_power_min_load
                * self.time_constant_onoff_operation
                * load_ratio
                * (1. - load_ratio)
                / time_constant_for_service;

            // Based on BS EN 15316-4-2:2017 equations 25 and 26, with
            // addition of weighting factor to handle existence of more than
            // one space heating service
            let compressor_power_part_load_average = compressor_power_full_load * load_ratio;

            // Based on BS EN 15316-4-2:2017 equation 29
            // Note: The variable t_H in equation 29 of EN 15316-4-2:2017 is
            // implied to be the running time required at full load, but
            // this is not consistent with the use of part-load compressor
            // power (P_gen,comp,LR) and leads to very low cycling energy
            // when multiplied by the cycling power (P_gen,comp,ONOF,LR) for
            // oversized heat pumps because the running time required at
            // full load for oversized heat pumps is relatively low.
            // Using the part-load running time (i.e. the time the heat pump
            // would run at the relevant load ratio) also does not make
            // sense given that both P_(gen,comp,LR) (from equations 25 and
            // 26)  and P_(gen,comp,ONOF,LR) (from equation 28) already
            // account for load ratio based on the time available after
            // water heating. In other words, both power variables are
            // already averaged over the time available rather than being
            // the power output only during the “on” part of the on-off
            // cycle.
            // Therefore, the approach used here is to use the time
            // available after running higher-priority services in equation
            // 29.
            let mut energy_input_hp_temp = (compressor_power_part_load_average
                * (1. + HEAT_PUMP_F_AUX)
                + power_used_due_to_inertia_effects)
                * time_available_for_current_service;

            // As the time available for space heating applies to all space
            // heating services combined, we need to avoid double-counting
            // by apportioning it between the space heating services in
            // proportion to the energy delivered for each space heating
            // service
            let weighting_factor = energy_delivered_hp / energy_delivered_hp_aggregated;

            energy_input_hp_temp *= weighting_factor;

            energy_input_hp_temp
        } else {
            // If not operating in on/off mode
            energy_delivered_hp / cop_op_cond
        };

        (energy_input_hp, compressor_power_min_load)
    }

    /// Calculate energy required by heat pump to satisfy demand for the service indicated.
    ///
    /// Note: Call via a HeatPumpService object, not directly.
    pub(crate) fn demand_energy(
        &mut self,
        service_name: &str,
        service_type: &HeatingServiceType,
        energy_output_required: f64,
        temp_output: Option<f64>,
        temp_return_feed: f64,
        temp_limit_upper: f64,
        design_flow_temp_op_cond: f64,
        time_constant_for_service: f64,
        service_on: bool,
        simtime: SimulationTimeIteration,
        temp_spread_correction: Option<TempSpreadCorrectionArg>,
        temp_used_for_scaling: Option<f64>,
        time_start: Option<f64>,
        hybrid_boiler_service: Option<HybridBoilerService>,
        emitters_data_for_buffer_tank: Option<BufferTankEmittersDataWithResult>,
        update_heat_source_state: Option<bool>,
    ) -> anyhow::Result<f64> {
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);
        let mut service_results = self.run_demand_energy_calc(
            service_name,
            service_type,
            energy_output_required,
            temp_output,
            temp_return_feed,
            temp_limit_upper,
            design_flow_temp_op_cond,
            time_constant_for_service,
            service_on,
            simtime,
            temp_spread_correction,
            temp_used_for_scaling,
            hybrid_boiler_service.clone(),
            None,
            None,
            Some(time_start.unwrap_or(0.)),
            emitters_data_for_buffer_tank,
        )?;
        let time_running = service_results.time_running_full_load;
        let energy_delivered_total = service_results.energy_delivered_total;

        let energy_output_delivered_boiler = if let Some(hybrid_boiler) = hybrid_boiler_service {
            // Call demand function for boiler and
            // return the boiler running time so that it can be added to the total running time
            let is_hybrid_service = true;

            let time_elapsed_hp = match self.backup_ctrl {
                HeatPumpBackupControlType::TopUp => None,
                HeatPumpBackupControlType::Substitute | HeatPumpBackupControlType::None => {
                    Some(self.total_time_running_current_timestep_full_load)
                }
            };

            let (energy_output_delivered_boiler, time_running_boiler) = match hybrid_boiler {
                HybridBoilerService::Regular(boiler_regular) => {
                    boiler_regular.lock().demand_energy(
                        service_results.energy_output_required_boiler,
                        kelvin_to_celsius(
                            temp_output.ok_or_else(|| anyhow!("Expected temp_output to be set"))?,
                        )?,
                        Some(kelvin_to_celsius(temp_return_feed)?),
                        Some(is_hybrid_service),
                        time_elapsed_hp,
                        Some(update_heat_source_state),
                        simtime,
                    )?
                }
                HybridBoilerService::Space(boiler_space) => boiler_space.lock().demand_energy(
                    service_results.energy_output_required_boiler,
                    kelvin_to_celsius(
                        temp_output.ok_or_else(|| anyhow!("Expected temp_output to be set"))?,
                    )?,
                    Some(kelvin_to_celsius(temp_return_feed)?),
                    time_start,
                    Some(is_hybrid_service),
                    time_elapsed_hp,
                    Some(update_heat_source_state),
                    simtime,
                )?,
            };
            if self.backup_ctrl == HeatPumpBackupControlType::Substitute && update_heat_source_state
            {
                self.total_time_running_current_timestep_full_load += time_running_boiler
                    .expect("hybrid_service_bool was true so value expected here");
            }

            energy_output_delivered_boiler
        } else {
            0.0
        };

        service_results.energy_output_delivered_boiler = Some(energy_output_delivered_boiler);

        // Save results that are needed later (in the timestep_end function)
        if update_heat_source_state {
            self.service_results
                .write()
                .push(ServiceResult::Full(Box::new(service_results)));
            self.total_time_running_current_timestep_full_load += time_running;
        }

        Ok(energy_delivered_total + energy_output_delivered_boiler)
    }

    /// Return the cumulative running time and throughput factor (exhaust air HPs only)
    pub(crate) fn running_time_throughput_factor(
        &mut self,
        space_heat_running_time_cumulative: f64,
        service_name: &str,
        service_type: &HeatingServiceType,
        energy_output_required: f64,
        temp_output: f64,
        temp_return_feed: f64,
        temp_limit_upper: f64,
        design_flow_temp_op_cond: f64,
        time_constant_for_service: f64,
        service_on: bool,
        volume_heated_by_service: f64,
        simtime: SimulationTimeIteration,
        temp_spread_correction: Option<TempSpreadCorrectionArg>,
        time_start: Option<f64>,
    ) -> anyhow::Result<(f64, f64)> {
        let time_start = time_start.unwrap_or(0.);
        let service_results = self.run_demand_energy_calc(
            service_name,
            service_type,
            energy_output_required,
            Some(temp_output),
            temp_return_feed,
            temp_limit_upper,
            design_flow_temp_op_cond,
            time_constant_for_service,
            service_on,
            simtime,
            temp_spread_correction,
            None,
            None,
            None,
            Some(space_heat_running_time_cumulative),
            Some(time_start),
            None,
        )?;

        let throughput_factor = self.calc_throughput_factor(service_results.time_running_full_load);

        // Adjust throughput factor to match simplifying assumption that all
        // extra ventilation required due to space heating demand in the current
        // zone is assigned to current zone
        let throughput_factor_zone = (throughput_factor - 1.0)
            * self
            .volume_heated_all_services
            .expect("Volume heated all services was expected to be set. This is only set for heat pumps with exhaust air.")
            / volume_heated_by_service
            + 1.0;
        Ok((
            service_results.time_running_full_load,
            throughput_factor_zone,
        ))
    }

    fn calc_throughput_factor(&self, time_running: f64) -> f64 {
        let timestep = self.simulation_timestep;

        // Apply overventilation ratio to part of timestep where HP is running
        // to calculate throughput_factor.
        ((timestep - time_running) + self.overvent_ratio * time_running) / timestep
    }

    pub fn throughput_factor(&self) -> f64 {
        self.calc_throughput_factor(self.total_time_running_current_timestep_full_load)
    }

    fn calc_energy_input(&self, t_idx: usize) -> anyhow::Result<()> {
        // Aggregate space heating services
        // TODO (from Python) This is only necessary because the model cannot handle an
        //  emitter circuit that serves more than one zone. If/when this
        //  capability is added, there will no longer be separate space
        //  heating services for each zone and this aggregation can be
        //  removed as it will not be necessary. At that point, the other
        //  contents of this function could also be moved back to their
        //  original locations
        let space_service_aggregates = self
            .service_results
            .read()
            .iter()
            .filter_map(|result| match result {
                ServiceResult::Full(x) if x.service_type == HeatingServiceType::Space => Some(x),
                _ => None,
            })
            .fold(None, |acc, x| match acc {
                None => Some((
                    x.time_running_full_load,
                    x.time_available_for_load_ratio_calc,
                    x.energy_delivered_hp,
                )),
                Some((run_sum, avail_max, energy_sum)) => Some((
                    run_sum + x.time_running_full_load,
                    avail_max.max(x.time_available_for_load_ratio_calc),
                    energy_sum + x.energy_delivered_hp,
                )),
            });

        for service_data in self.service_results.write().iter_mut() {
            if let ServiceResult::Full(service_data) = service_data {
                let temp_output = service_data
                    .temp_output
                    .ok_or_else(|| anyhow!("Expected temp_output to be set"))?;
                let energy_input_backup = service_data.energy_input_backup;
                let energy_heating_circ_pump = service_data.energy_heating_circ_pump;
                let energy_source_circ_pump = service_data.energy_source_circ_pump;
                let energy_heating_warm_air_fan = service_data.energy_heating_warm_air_fan;
                let use_backup_heater_only = service_data.use_backup_heater_only;

                let (
                    load_ratio,
                    load_ratio_continuous_min,
                    hp_operating_in_onoff_mode,
                    time_available_for_current_service,
                    energy_delivered_hp_aggregated,
                ) = match service_data.service_type {
                    HeatingServiceType::DomesticHotWaterRegular => (
                        1.,
                        1.,
                        false,
                        service_data.time_available_for_load_ratio_calc,
                        service_data.energy_delivered_hp,
                    ),
                    HeatingServiceType::Space => {
                        let (
                            time_running_for_load_ratio,
                            time_available_for_current_service,
                            energy_delivered_hp_aggregated,
                        ) = space_service_aggregates
                            .expect("Space service found but aggregates missing");
                        // TODO (from Python) Check that certain parameters are the same across all space heating services
                        let (load_ratio, load_ratio_continuous_min, hp_operating_in_onoff_mode) =
                            self.load_ratio_and_mode(
                                time_running_for_load_ratio,
                                time_available_for_current_service,
                                temp_output,
                            )?;

                        (
                            load_ratio,
                            load_ratio_continuous_min,
                            hp_operating_in_onoff_mode,
                            time_available_for_current_service,
                            energy_delivered_hp_aggregated,
                        )
                    }
                    _ => bail!("Service type should be either DomesticHotWaterRegular or Space"),
                };

                if use_backup_heater_only {
                    service_data.time_running_part_load = Some(service_data.time_running_full_load)
                } else {
                    service_data.time_running_part_load = Some(
                        service_data.time_running_full_load
                            / max_of_2(load_ratio, load_ratio_continuous_min),
                    );
                }

                let (energy_input_hp, compressor_power_min_load) = self.energy_input_compressor(
                    service_data.service_on,
                    service_data.use_backup_heater_only,
                    hp_operating_in_onoff_mode,
                    service_data.energy_delivered_hp,
                    energy_delivered_hp_aggregated,
                    service_data.thermal_capacity_op_cond,
                    service_data
                        .cop_op_cond
                        .ok_or_else(|| anyhow::anyhow!("Expected cop_op_cond to be set"))?,
                    time_available_for_current_service,
                    load_ratio,
                    load_ratio_continuous_min,
                    service_data.time_constant_for_service,
                    service_data.service_type,
                );

                let energy_input_total = energy_input_hp
                    + energy_input_backup
                    + energy_heating_circ_pump
                    + energy_source_circ_pump
                    + energy_heating_warm_air_fan;

                service_data.compressor_power_min_load = compressor_power_min_load;
                service_data.load_ratio_continuous_min = load_ratio_continuous_min;
                service_data.load_ratio = load_ratio;
                service_data.hp_operating_in_onoff_mode = hp_operating_in_onoff_mode;
                service_data.energy_input_hp = energy_input_hp;
                service_data.energy_input_total = energy_input_total;

                // Feed/return results to other modules
                self.energy_supply_connections[&service_data.service_name]
                    .demand_energy(energy_input_total, t_idx)?;
            }
        }

        Ok(())
    }

    /// Calculate auxiliary energy according to CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.7
    fn calc_auxiliary_energy(
        &self,
        timestep: f64,
        time_remaining_current_timestep_part_load: f64,
        timestep_idx: usize,
    ) -> anyhow::Result<(f64, f64, f64)> {
        // Retrieve control settings for this timestep
        let mut heating_profile_on = false;
        let mut water_profile_on = false;
        for service_data in self.service_results.read().iter() {
            if let ServiceResult::Full(calculation) = service_data {
                match calculation.service_type {
                    HeatingServiceType::Space => {
                        heating_profile_on = calculation.service_on;
                    }
                    HeatingServiceType::DomesticHotWaterRegular => {
                        water_profile_on = calculation.service_on;
                    }
                    _ => bail!("Service type should be either DomesticHotWaterRegular or Space"),
                }
            }
        }

        // Energy used in standby and crankcase heater mode
        let mut energy_off_mode = 0.0;
        let mut energy_standby = 0.0;
        let mut energy_crankcase_heater_mode = 0.0;
        match (heating_profile_on, water_profile_on) {
            (true, _) => {
                energy_standby = time_remaining_current_timestep_part_load * self.power_standby;
                energy_crankcase_heater_mode =
                    time_remaining_current_timestep_part_load * self.power_crankcase_heater_mode;
            }
            (false, true) => {
                energy_standby = time_remaining_current_timestep_part_load * self.power_standby;
            }
            (false, false) => {
                energy_off_mode = timestep * self.power_off_mode;
            }
        }

        let energy_aux = energy_standby + energy_crankcase_heater_mode + energy_off_mode;
        self.energy_supply_connection_aux
            .demand_energy(energy_aux, timestep_idx)?;

        Ok((
            energy_standby,
            energy_crankcase_heater_mode,
            energy_off_mode,
        ))
    }

    /// If HP uses heat network as source, calculate energy extracted from heat network
    fn extract_energy_from_source(&self, timestep_idx: usize) -> anyhow::Result<()> {
        for service_data in self.service_results.read().iter() {
            if let ServiceResult::Full(service_data) = service_data {
                let HeatPumpEnergyCalculation {
                    service_name,
                    energy_delivered_hp,
                    energy_input_hp,
                    ..
                } = service_data.as_ref();
                let energy_extracted_hp = energy_delivered_hp - energy_input_hp;

                if energy_extracted_hp < 0.
                    && !is_close!(energy_extracted_hp, 0., rel_tol = 1e-09, abs_tol = 1e-10)
                {
                    bail!("Energy extracted from source ({energy_extracted_hp}) by heat pump should not be negative for service: {service_name}")
                }
                self.energy_supply_heat_source_connections[service_name.as_str()]
                    .demand_energy(energy_extracted_hp, timestep_idx)?;
            }
        }
        Ok(())
    }

    /// Calculations to be done at the end of each timestep
    pub fn timestep_end(&mut self, timestep_idx: usize) -> anyhow::Result<()> {
        self.calc_energy_input(timestep_idx)?;

        let timestep = self.simulation_timestep;
        let time_remaining_current_timestep_full_load =
            timestep - self.total_time_running_current_timestep_full_load;
        let part_load_sum: f64 = self
            .service_results
            .read()
            .iter()
            .filter_map(|result| {
                if let ServiceResult::Full(heat_pump_calc) = result {
                    heat_pump_calc.time_running_part_load
                } else {
                    None
                }
            })
            .sum();
        let time_remaining_current_timestep_part_load = timestep - part_load_sum;

        if time_remaining_current_timestep_full_load == 0.0 {
            self.time_running_continuous += self.total_time_running_current_timestep_full_load;
        } else {
            self.time_running_continuous = 0.;
        }

        let (energy_standby, energy_crankcase_heater_mode, energy_off_mode) = self
            .calc_auxiliary_energy(
                timestep,
                time_remaining_current_timestep_part_load,
                timestep_idx,
            )?;

        if self.energy_supply_heat_source.is_some() {
            self.extract_energy_from_source(timestep_idx)?;
        }

        // If detailed results are to be output, save the results from the current timestep
        if let Some(ref mut detailed_results) = self.detailed_results {
            self.service_results
                .write()
                .push(ServiceResult::Aux(AuxiliaryParameters {
                    energy_standby,
                    energy_crankcase_heater_mode,
                    energy_off_mode,
                }));
            detailed_results.push(self.service_results.clone());
        }

        // Variables below need to be reset at the end of each timestep.
        self.total_time_running_current_timestep_full_load = Default::default();
        self.service_results = Default::default();

        Ok(())
    }

    pub fn output_detailed_results(
        &self,
        hot_water_energy_output: &IndexMap<String, Vec<ResultParamValue>>,
        hot_water_source_name_for_heat_pump_service: &IndexMap<String, String>,
    ) -> (ResultsPerTimestep, ResultsAnnual) {
        let detailed_results = self.detailed_results.as_ref().expect(
            "Detailed results cannot be output when the option to collect them was not selected",
        );
        let mut results_per_timestep: ResultsPerTimestep = Default::default();
        results_per_timestep.insert("auxiliary".into(), Default::default());
        for (parameter, param_unit, _) in AUX_PARAMETERS {
            let auxiliary_param_results = results_per_timestep
                .get_mut("auxiliary")
                .unwrap()
                .entry((parameter.into(), Some(param_unit.into())))
                .or_default();
            for service_results in detailed_results {
                if let ServiceResult::Aux(calc) = service_results
                    .read()
                    .iter()
                    .last()
                    .expect("Expected at least one service result")
                {
                    let result = calc.param(parameter);
                    auxiliary_param_results.push(result);
                }
            }
        }

        let energy_supply_connections = self.energy_supply_connections.clone();

        // For each service, report required output parameters
        for (service_idx, service_name) in energy_supply_connections.keys().enumerate() {
            let results_entry = results_per_timestep
                .entry(service_name.clone())
                .or_default();
            // Look up each required parameter
            for (parameter, param_unit, _) in OUTPUT_PARAMETERS {
                let param_entry = results_entry
                    .entry((parameter.into(), param_unit.map(Into::into)))
                    .or_default();
                // Look up value of required parameter in each timestep
                for service_results in self
                    .detailed_results
                    .as_ref()
                    .expect("Detailed results are expected to have been initialised")
                    .iter()
                {
                    let result = match &service_results.read()[service_idx] {
                        ServiceResult::Full(calc) => calc.param(parameter),
                        ServiceResult::Aux(_) => unreachable!(),
                    };
                    param_entry.push(result);
                }
            }
            // For water heating service, record hot water energy delivered from tank
            if matches!(
                match &self.detailed_results.as_ref().unwrap()[0].read()[service_idx] {
                    ServiceResult::Full(calc) => calc.service_type,
                    ServiceResult::Aux(_) => unreachable!(),
                },
                HeatingServiceType::DomesticHotWaterRegular
            ) {
                let energy_delivered_total_len = results_per_timestep[service_name]
                    [&("energy_delivered_total".into(), Some("kWh".into()))]
                    .len();
                results_per_timestep.get_mut(service_name).unwrap().insert(
                    ("energy_delivered_H5".into(), Some("kWh".into())),
                    {
                        // For DHW, need to include storage and primary circuit losses.
                        // Can do this by replacing H5 numerator with total energy
                        // draw-off from hot water cylinder.
                        let hws_name = &hot_water_source_name_for_heat_pump_service[service_name];
                        if !hot_water_energy_output.contains_key(hws_name) {
                            vec![ResultParamValue::Empty; energy_delivered_total_len]
                        } else {
                            hot_water_energy_output[hws_name].clone()
                        }
                    },
                );
            } else {
                let energy_delivered_total_results = results_per_timestep[service_name]
                    [&("energy_delivered_total".into(), Some("kWh".into()))]
                    .clone();
                results_per_timestep.get_mut(service_name).unwrap().insert(
                    ("energy_delivered_H5".into(), Some("kWh".into())),
                    energy_delivered_total_results,
                );
            }
        }

        let mut results_annual: ResultsAnnual = Default::default();
        results_annual.insert(
            "Overall".into(),
            OUTPUT_PARAMETERS
                .iter()
                .filter(|param| param.2)
                .map(|(parameter, param_units, _)| {
                    (
                        ((*parameter).into(), param_units.map(Into::into).to_owned()),
                        ResultParamValue::Number(0.0),
                    )
                })
                .collect::<ResultAnnual>(),
        );
        results_annual.insert("auxiliary".into(), Default::default());
        results_annual.get_mut("Overall").unwrap().insert(
            ("energy_delivered_H5".into(), Some("kWh".into())),
            ResultParamValue::Number(0.0),
        );
        // Report auxiliary parameters (not specific to a service)
        {
            // make scope for this mutable borrow
            let auxiliary_annual_results = results_annual.get_mut("auxiliary").unwrap();
            for (parameter, param_unit, incl_in_annual) in AUX_PARAMETERS {
                if incl_in_annual {
                    auxiliary_annual_results.insert(
                        (parameter.into(), Some(param_unit.into())),
                        results_per_timestep["auxiliary"]
                            [&(parameter.into(), Some(param_unit.into()))]
                            .iter()
                            .cloned()
                            .sum::<ResultParamValue>(),
                    );
                }
            }
        }

        // For each service, report required output parameters
        for service_name in energy_supply_connections.keys() {
            let param_totals_for_overall = {
                let annual_results_entry = results_annual.entry(service_name.clone()).or_default();
                let mut param_totals_for_overall: HashMap<(&str, &str), ResultParamValue> =
                    Default::default();
                for (parameter, param_unit, incl_in_annual) in OUTPUT_PARAMETERS {
                    if incl_in_annual {
                        let parameter_annual_total = results_per_timestep[service_name.as_str()]
                            [&(parameter.into(), param_unit.map(Into::into))]
                            .iter()
                            .cloned()
                            .sum::<ResultParamValue>();
                        annual_results_entry.insert(
                            (parameter.into(), param_unit.map(Into::into)),
                            parameter_annual_total.clone(),
                        );
                        param_totals_for_overall.insert(
                            (parameter, param_unit.unwrap_or("")),
                            parameter_annual_total,
                        );
                    }
                }

                param_totals_for_overall
            };
            for (param_index, param_total_for_overall) in param_totals_for_overall {
                *results_annual
                    .get_mut("Overall")
                    .unwrap()
                    .get_mut(&(param_index.0.into(), Some(param_index.1.into())))
                    .unwrap() += param_total_for_overall;
            }

            if results_per_timestep[service_name]
                [&("energy_delivered_H5".into(), Some("kWh".into()))]
                .contains(&ResultParamValue::Empty)
            {
                results_annual.get_mut(service_name).unwrap().insert(
                    ("energy_delivered_H5".into(), Some("kWh".into())),
                    ResultParamValue::Empty,
                );
                *results_annual
                    .get_mut("Overall")
                    .unwrap()
                    .get_mut(&("energy_delivered_H5".into(), Some("kWh".into())))
                    .unwrap() = ResultParamValue::Empty;
            } else {
                results_annual.get_mut(service_name).unwrap().insert(
                    ("energy_delivered_H5".into(), Some("kWh".into())),
                    results_per_timestep[service_name.as_str()]
                        [&("energy_delivered_H5".into(), Some("kWh".into()))]
                        .iter()
                        .cloned()
                        .sum::<ResultParamValue>(),
                );

                if results_annual["Overall"][&("energy_delivered_H5".into(), Some("kWh".into()))]
                    != ResultParamValue::Empty
                {
                    let service_energy_delivered_results = results_annual[service_name.as_str()]
                        [&("energy_delivered_H5".into(), Some("kWh".into()))]
                        .clone();

                    *results_annual
                        .get_mut("Overall")
                        .unwrap()
                        .get_mut(&("energy_delivered_H5".into(), Some("kWh".into())))
                        .unwrap() += service_energy_delivered_results;
                }
            }

            self.calc_service_cop(results_annual.get_mut(service_name).unwrap(), None);
        }

        // Calculate overall CoP for all services combined
        self.calc_service_cop_overall(&mut results_annual);

        (results_per_timestep, results_annual)
    }

    /// Calculate CoP for whole simulation period for the given service (or overall)
    fn calc_service_cop_overall(&self, results_annual: &mut ResultsAnnual) {
        let results_auxiliary = results_annual.get("auxiliary").cloned();
        let results_totals = results_annual.get_mut("Overall").unwrap();

        self.calc_service_cop(results_totals, results_auxiliary);
    }

    fn calc_service_cop(
        &self,
        results_totals: &mut ResultAnnual,
        results_auxiliary: Option<ResultAnnual>,
    ) {
        // Add auxiliary energy to overall CoP
        let energy_auxiliary = {
            match results_auxiliary {
                Some(results_aux) => results_aux.values().cloned().sum(),
                None => ResultParamValue::Number(0.),
            }
        };

        // Calculate CoP at different system boundaries
        let cop_h1_numerator =
            results_totals[&("energy_delivered_HP".into(), Some("kWh".into()))].clone();
        let cop_h1_denominator =
            &results_totals[&("energy_input_HP".into(), Some("kWh".into()))] + &energy_auxiliary;
        let cop_h2_numerator = cop_h1_numerator.clone();
        let cop_h2_denominator = &cop_h1_denominator
            + &results_totals[&("energy_source_circ_pump".into(), Some("kWh".into()))];
        let cop_h3_numerator = &cop_h2_numerator
            + &results_totals[&("energy_delivered_backup".into(), Some("kWh".into()))].clone();
        let cop_h3_denominator = &cop_h2_denominator
            + &results_totals[&("energy_input_backup".into(), Some("kWh".into()))].clone();
        let cop_h4_numerator = cop_h3_numerator.clone();
        let cop_h4_denominator = &cop_h3_denominator
            + &results_totals[&("energy_heating_circ_pump".into(), Some("kWh".into()))].clone();
        let cop_h5_numerator =
            results_totals[&("energy_delivered_H5".into(), Some("kWh".into()))].clone();
        let cop_h5_denominator = cop_h4_denominator.clone();

        results_totals.insert(
            ("CoP (H1)".into(), None),
            if cop_h1_denominator == 0.0 {
                0.0.into()
            } else {
                cop_h1_numerator / cop_h1_denominator
            },
        );
        results_totals.insert(
            ("CoP (H2)".into(), None),
            if cop_h2_denominator == 0.0 {
                0.0.into()
            } else {
                cop_h2_numerator / cop_h2_denominator
            },
        );
        results_totals.insert(
            ("CoP (H3)".into(), None),
            if cop_h3_denominator == 0.0 {
                0.0.into()
            } else {
                cop_h3_numerator / cop_h3_denominator
            },
        );
        results_totals.insert(
            ("CoP (H4)".into(), None),
            if cop_h4_denominator == 0.0 {
                0.0.into()
            } else {
                cop_h4_numerator / cop_h4_denominator
            },
        );

        let (subkey, value) = if cop_h5_denominator == 0. {
            (None, ResultParamValue::from(0.))
        } else if cop_h5_numerator == ResultParamValue::Empty {
            let cop_h5_note = "Note: Cannot calculate CoP (H5) when HP is heating a pre-heat tank";

            (Some(cop_h5_note.into()), ResultParamValue::Empty)
        } else {
            let cop_h5_note =
                "Note: For water heating services, only valid when HP is only heat source";

            (
                Some(cop_h5_note.into()),
                cop_h5_numerator / cop_h5_denominator,
            )
        };
        results_totals.insert(("CoP (H5)".into(), subkey), value);
    }
}

const OUTPUT_PARAMETERS: [(&str, Option<&str>, bool); 22] = [
    ("service_name", None, false),
    ("service_type", None, false),
    ("service_on", None, false),
    ("energy_output_required", Some("kWh"), true),
    ("temp_output", Some("K"), false),
    ("temp_source", Some("K"), false),
    ("thermal_capacity_op_cond", Some("kW"), false),
    ("cop_op_cond", None, false),
    ("time_running_full_load", Some("hours"), true),
    ("time_running_part_load", Some("hours"), true),
    ("load_ratio", None, false),
    ("hp_operating_in_onoff_mode", None, false),
    ("energy_delivered_HP", Some("kWh"), true),
    ("energy_delivered_backup", Some("kWh"), true),
    ("energy_delivered_total", Some("kWh"), true),
    ("energy_input_HP", Some("kWh"), true),
    ("energy_input_backup", Some("kWh"), true),
    ("energy_heating_circ_pump", Some("kWh"), true),
    ("energy_source_circ_pump", Some("kWh"), true),
    ("energy_heating_warm_air_fan", Some("kWh"), true),
    ("energy_input_total", Some("kWh"), true),
    ("energy_output_delivered_boiler", Some("kWh"), true),
];
const AUX_PARAMETERS: [(&str, &str, bool); 3] = [
    ("energy_standby", "kWh", true),
    ("energy_crankcase_heater_mode", "kWh", true),
    ("energy_off_mode", "kWh", true),
];

type ResultAnnual = IndexMap<(String, Option<String>), ResultParamValue>;

fn result_str(string: &str) -> String {
    string.into()
}

pub(crate) enum TempSpreadCorrectionArg {
    Float(f64),
    Callable(Box<dyn Fn(f64, f64, f64, HeatPumpTestData) -> f64>),
}

impl Debug for TempSpreadCorrectionArg {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                TempSpreadCorrectionArg::Float(f) => format!("float: {f}"),
                TempSpreadCorrectionArg::Callable(_) =>
                    "A callable taking three f64 values and returning an f64.".to_owned(),
            }
        )
    }
}

#[derive(Clone)]
pub(crate) enum HybridBoilerService {
    Regular(Arc<Mutex<BoilerServiceWaterRegular>>),
    Space(Arc<Mutex<BoilerServiceSpace>>),
}

#[derive(Clone, Debug, PartialEq)]
pub enum ServiceResult {
    Full(Box<HeatPumpEnergyCalculation>),
    Aux(AuxiliaryParameters),
}

#[derive(Clone, Debug, PartialEq)]
pub struct HeatPumpEnergyCalculation {
    service_name: String,
    service_type: HeatingServiceType,
    service_on: bool,
    energy_output_required: f64,
    temp_output: Option<f64>,
    temp_source: f64,
    cop_op_cond: Option<f64>,
    thermal_capacity_op_cond: Option<f64>,
    time_running_full_load: f64,
    time_running_part_load: Option<f64>,
    time_available_for_capacity_calc: f64,
    time_available_for_load_ratio_calc: f64,
    time_constant_for_service: f64,
    compressor_power_min_load: f64,
    load_ratio_continuous_min: f64,
    load_ratio: f64,
    use_backup_heater_only: bool,
    hp_operating_in_onoff_mode: bool,
    energy_input_hp: f64,
    energy_delivered_hp: f64,
    energy_input_backup: f64,
    energy_delivered_backup: f64,
    energy_input_total: f64,
    energy_delivered_total: f64,
    energy_heating_circ_pump: f64,
    energy_source_circ_pump: f64,
    energy_output_required_boiler: f64,
    energy_heating_warm_air_fan: f64,
    energy_output_delivered_boiler: Option<f64>, // this field is usually set later - not splitting out new type for now
}

impl HeatPumpEnergyCalculation {
    fn param(&self, param: &str) -> ResultParamValue {
        match param {
            "service_name" => ResultParamValue::String(self.service_name.clone()),
            "service_type" => ResultParamValue::String(self.service_type.to_string().into()),
            "service_on" => ResultParamValue::Boolean(self.service_on),
            "energy_output_required" => ResultParamValue::Number(self.energy_output_required),
            "temp_output" => ResultParamValue::Number(
                self.temp_output
                    .expect("The temp output value is expected to have been calculated."),
            ),
            "temp_source" => ResultParamValue::Number(self.temp_source),
            "thermal_capacity_op_cond" => ResultParamValue::Number(
                self.thermal_capacity_op_cond
                    .expect("The thermal capacity op cond is expected to have been calculated."),
            ),
            "cop_op_cond" => ResultParamValue::Number(
                self.cop_op_cond
                    .expect("The cop op cond value is expected to have been calculated."),
            ),
            "time_running_full_load" => ResultParamValue::Number(self.time_running_full_load),
            "time_running_part_load" => {
                ResultParamValue::Number(self.time_running_part_load.expect(
                    "The time_running_part_load value is expected to have been calculated.",
                ))
            }
            "load_ratio" => ResultParamValue::Number(self.load_ratio),
            "hp_operating_in_onoff_mode" => {
                ResultParamValue::Boolean(self.hp_operating_in_onoff_mode)
            }
            "energy_delivered_HP" => ResultParamValue::Number(self.energy_delivered_hp),
            "energy_delivered_backup" => ResultParamValue::Number(self.energy_delivered_backup),
            "energy_delivered_total" => ResultParamValue::Number(self.energy_delivered_total),
            "energy_input_HP" => ResultParamValue::Number(self.energy_input_hp),
            "energy_input_backup" => ResultParamValue::Number(self.energy_input_backup),
            "energy_heating_circ_pump" => ResultParamValue::Number(self.energy_heating_circ_pump),
            "energy_source_circ_pump" => ResultParamValue::Number(self.energy_source_circ_pump),
            "energy_heating_warm_air_fan" => {
                ResultParamValue::Number(self.energy_heating_warm_air_fan)
            }
            "energy_input_total" => ResultParamValue::Number(self.energy_input_total),
            "energy_output_delivered_boiler" => {
                ResultParamValue::Number(self.energy_output_delivered_boiler.expect(
                    "The energy output delivered boiler value is expected to have been calculated.",
                ))
            }
            &_ => panic!("Parameter {param} not recognised"),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AuxiliaryParameters {
    energy_standby: f64,
    energy_crankcase_heater_mode: f64,
    energy_off_mode: f64,
}

impl AuxiliaryParameters {
    pub(crate) fn param(&self, param: &str) -> ResultParamValue {
        match param {
            "energy_standby" => ResultParamValue::Number(self.energy_standby),
            "energy_crankcase_heater_mode" => {
                ResultParamValue::Number(self.energy_crankcase_heater_mode)
            }
            "energy_off_mode" => ResultParamValue::Number(self.energy_off_mode),
            &_ => panic!("Parameter {param} not recognised"),
        }
    }
}

#[derive(Clone, Debug)]
struct CostScheduleMetadata {
    pub hp: Vec<f64>,
    pub boiler: Vec<f64>,
    pub start_day: u32,
    pub time_series_step: f64,
}

impl TryFrom<&BoilerCostScheduleHybrid> for CostScheduleMetadata {
    type Error = anyhow::Error;

    fn try_from(cost_schedule: &BoilerCostScheduleHybrid) -> anyhow::Result<Self> {
        Ok(Self {
            hp: reject_nulls(expand_numeric_schedule(&cost_schedule.cost_schedule_hp))?,
            boiler: reject_nulls(expand_numeric_schedule(&cost_schedule.cost_schedule_boiler))?,
            start_day: cost_schedule.cost_schedule_start_day,
            time_series_step: cost_schedule.cost_schedule_time_series_step,
        })
    }
}

// In Python this is HeatPump_HWOnly
/// An object to represent an electric hot-water-only heat pump, tested to EN 16147
#[derive(Clone, Debug)]
pub struct HeatPumpHotWaterOnly {
    power_in_kw: f64,
    energy_supply_connection: EnergySupplyConnection,
    simulation_timestep: f64,
    control_min: Arc<Control>,
    control_max: Arc<Control>,
    initial_efficiency: f64,
    tank_volume: f64,
    daily_losses: f64,
    heat_exchanger_surface_area: f64,
    in_use_factor_mismatch: f64,
    tank_volume_declared: f64,
    heat_exchanger_surface_area_declared: f64,
    daily_losses_declared: f64,
}

impl HeatPumpHotWaterOnly {
    pub(crate) fn new(
        power_max_in_kw: f64,
        energy_supply_connection: EnergySupplyConnection,
        test_data: &HeatPumpHotWaterTestData,
        vol_daily_average: f64,
        tank_volume: f64,
        daily_losses: f64,
        heat_exchanger_surface_area: f64,
        in_use_factor_mismatch: f64,
        tank_volume_declared: f64,
        heat_exchanger_surface_area_declared: f64,
        daily_losses_declared: f64,
        simulation_timestep: f64,
        control_min: Arc<Control>, // TODO in Python 1.0.0a1 this is TimeControl
        control_max: Arc<Control>, // TODO in Python 1.0.0a1 this is TimeControl
    ) -> Self {
        let efficiencies = Efficiencies {
            l: test_data.l.as_ref().map(|profile_data| {
                let HeatPumpHotWaterOnlyTestDatum {
                    cop_dhw,
                    hw_tapping_prof_daily_total: hw_tapping_prof_daily,
                    energy_input_measured,
                    power_standby,
                    hw_vessel_loss_daily,
                } = *profile_data;
                Self::init_efficiency_tapping_profile(
                    cop_dhw,
                    hw_tapping_prof_daily,
                    energy_input_measured,
                    power_standby,
                    hw_vessel_loss_daily,
                )
            }),
            m: {
                let HeatPumpHotWaterOnlyTestDatum {
                    cop_dhw,
                    hw_tapping_prof_daily_total: hw_tapping_prof_daily,
                    energy_input_measured,
                    power_standby,
                    hw_vessel_loss_daily,
                } = test_data.m;
                Self::init_efficiency_tapping_profile(
                    cop_dhw,
                    hw_tapping_prof_daily,
                    energy_input_measured,
                    power_standby,
                    hw_vessel_loss_daily,
                )
            },
        };

        Self {
            power_in_kw: power_max_in_kw,
            energy_supply_connection,
            simulation_timestep,
            control_min,
            control_max,
            initial_efficiency: Self::init_efficiency(&efficiencies, vol_daily_average),
            tank_volume,
            daily_losses,
            heat_exchanger_surface_area,
            in_use_factor_mismatch,
            tank_volume_declared,
            heat_exchanger_surface_area_declared,
            daily_losses_declared,
        }
    }

    /// Calculate efficiency for given test condition (tapping profile)
    fn init_efficiency_tapping_profile(
        cop_dhw: f64,
        hw_tapping_prof_daily_total: f64,
        energy_input_measured: f64,
        power_standby: f64,
        hw_vessel_loss_daily: f64,
    ) -> f64 {
        // CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.2
        let temp_factor = 0.6 * 0.9;
        let energy_input_hw_vessel_loss = hw_vessel_loss_daily / cop_dhw * temp_factor;
        let energy_input_standby = power_standby * HOURS_PER_DAY as f64 * temp_factor;
        let energy_input_test =
            energy_input_measured - energy_input_standby + energy_input_hw_vessel_loss;
        let energy_demand_test = hw_tapping_prof_daily_total + hw_vessel_loss_daily * temp_factor;

        energy_demand_test / energy_input_test
    }

    /// Calculate initial efficiency based on SAP 10.2 section N3.7 b) and c)
    fn init_efficiency(efficiencies: &Efficiencies, vol_daily_average: f64) -> f64 {
        let eff_m = efficiencies.m;
        match efficiencies.l {
            None => efficiencies.m,
            Some(eff_l) => {
                let vol_daily_limit_lower = 100.2;
                let vol_daily_limit_upper = 199.8;
                if vol_daily_average <= vol_daily_limit_lower {
                    eff_m
                } else if vol_daily_average >= vol_daily_limit_upper {
                    // TODO report bug in Python in the conditional for this clause where
                    // self.__vol_daily_average is referenced without it having been init'd first
                    eff_l
                } else {
                    eff_m
                        + (eff_l - eff_m) / (vol_daily_limit_upper - vol_daily_limit_lower)
                            * (vol_daily_average - vol_daily_limit_lower)
                }
            }
        }
    }

    /// Calculate efficiency after applying in use factor if entered tank characteristics
    /// do not meet criteria of data in database.
    fn calc_efficiency(&self) -> f64 {
        let in_use_factor_mismatch = if self.tank_volume < self.tank_volume_declared
            || self.heat_exchanger_surface_area < self.heat_exchanger_surface_area_declared
            || self.daily_losses > self.daily_losses_declared
        {
            self.in_use_factor_mismatch
        } else {
            // no in use factor is applied
            1.0
        };

        self.initial_efficiency * in_use_factor_mismatch
    }

    pub(crate) fn setpnt(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, Option<f64>)> {
        // TODO review - this check may be able to be removed in future if we validate control earlier
        match (&self.control_min.as_ref(), &self.control_max.as_ref()) {
            (
                Control::CombinationTime { .. } | Control::SetpointTime { .. },
                Control::CombinationTime { .. } | Control::SetpointTime { .. },
            ) => Ok((
                self.control_min.setpnt(&simulation_time_iteration),
                self.control_max.setpnt(&simulation_time_iteration),
            )),
            _ => bail!(
                "Expected control_min and control_max to be combination or setpoint time controls"
            ),
        }
    }

    /// Demand energy (in kWh) from the heat pump
    pub fn demand_energy(
        &self,
        energy_demand: f64,
        _temp_flow: f64,
        _temp_return: f64,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        // Account for time control. In the Python they also check here whether control_min is None
        // but this is not possible as it's a required field for a HeatPumpHotWaterOnly object.
        let energy_supplied = if self.control_min.as_ref().is_on(simtime) {
            min_of_2(energy_demand, self.power_in_kw * self.simulation_timestep)
        } else {
            0.0
        };

        let energy_required = energy_supplied / self.calc_efficiency();
        self.energy_supply_connection
            .demand_energy(energy_required, simtime.index)
            .unwrap();

        energy_supplied
    }

    /// Calculate the maximum energy output (in kWh) from the heater
    pub fn energy_output_max(&self, _temp_return: f64, simtime: SimulationTimeIteration) -> f64 {
        // Account for time control. In the Python they also check here whether control_min is None
        // but this is not possible as it's a required field for a HeatPumpHotWaterOnly object.
        if self.control_min.as_ref().is_on(simtime) {
            self.power_in_kw * self.simulation_timestep
        } else {
            0.0
        }
    }
}

#[derive(Default)]
struct Efficiencies {
    l: Option<f64>,
    m: f64,
}

#[cfg(test)]
mod tests {
    use std::ops::Deref;

    use super::*;
    use crate::core::controls::time_control::{OnOffTimeControl, SetpointTimeControl};
    use crate::core::energy_supply::energy_supply::EnergySupplyBuilder;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::external_conditions::DaylightSavingsConfig;
    use crate::input::{BoilerHotWaterTest, FuelType, HeatPumpBufferTank, HeatSourceLocation};
    use crate::simulation_time::SimulationTime;
    use crate::{core::material_properties::WATER, external_conditions::ShadingSegment, input};
    use approx::{assert_relative_eq, assert_ulps_eq};
    use indexmap::indexmap;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use serde_json::{json, Value};

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[fixture]
    fn buffer_tank(simulation_time: SimulationTime) -> BufferTank {
        BufferTank::new(
            1.68,
            50.,
            15.,
            0.040,
            2,
            simulation_time.step,
            20.,
            *WATER,
            Some(true),
        )
    }

    // emitters data needed to run buffer tank calculations
    const EMITTERS_DATA_FOR_BUFFER_TANK: [BufferTankEmittersData; 8] = [
        BufferTankEmittersData {
            temp_emitter_req: 43.32561228292832,
            power_req_from_buffer_tank: 6.325422354229758,
            design_flow_temp: 55.,
            target_flow_temp: 48.54166666666667,
            temp_rm_prev: 22.488371468978006,
            variable_flow: true,
            temp_diff_emit_design: 10.,
            min_flow_rate: 0.05,
            max_flow_rate: 0.3,
        },
        BufferTankEmittersData {
            temp_emitter_req: 30.778566169260767,
            power_req_from_buffer_tank: 4.036801431329545,
            design_flow_temp: 55.,
            target_flow_temp: 48.54166666666667,
            temp_rm_prev: 22.61181775388348,
            variable_flow: true,
            temp_diff_emit_design: 10.,
            min_flow_rate: 0.05,
            max_flow_rate: 0.3,
        },
        BufferTankEmittersData {
            temp_emitter_req: 50.248361664005266,
            power_req_from_buffer_tank: 8.808940459526795,
            design_flow_temp: 55.,
            target_flow_temp: 49.375,
            temp_rm_prev: 17.736483875769345,
            variable_flow: true,
            temp_diff_emit_design: 10.,
            min_flow_rate: 0.05,
            max_flow_rate: 0.3,
        },
        BufferTankEmittersData {
            temp_emitter_req: 30.723032018863076,
            power_req_from_buffer_tank: 6.444185473658857,
            design_flow_temp: 55.,
            target_flow_temp: 49.375,
            temp_rm_prev: 16.85636993835381,
            variable_flow: true,
            temp_diff_emit_design: 10.,
            min_flow_rate: 0.05,
            max_flow_rate: 0.3,
        },
        BufferTankEmittersData {
            temp_emitter_req: 29.646908204800425,
            power_req_from_buffer_tank: 5.00550934831923,
            design_flow_temp: 55.,
            target_flow_temp: 48.75,
            temp_rm_prev: 17.22290169647781,
            variable_flow: true,
            temp_diff_emit_design: 10.,
            min_flow_rate: 0.05,
            max_flow_rate: 0.3,
        },
        BufferTankEmittersData {
            temp_emitter_req: 30.42417311801241,
            power_req_from_buffer_tank: 2.5923631086027457,
            design_flow_temp: 55.,
            target_flow_temp: 48.22916666666667,
            temp_rm_prev: 21.897823675853978,
            variable_flow: true,
            temp_diff_emit_design: 10.,
            min_flow_rate: 0.05,
            max_flow_rate: 0.3,
        },
        BufferTankEmittersData {
            temp_emitter_req: 39.38823309970535,
            power_req_from_buffer_tank: 1.9107311055974638,
            design_flow_temp: 55.,
            target_flow_temp: 48.4375,
            temp_rm_prev: 22.36470513987037,
            variable_flow: true,
            temp_diff_emit_design: 10.,
            min_flow_rate: 0.05,
            max_flow_rate: 0.3,
        },
        BufferTankEmittersData {
            temp_emitter_req: 28.934946810199524,
            power_req_from_buffer_tank: 1.711635686322881,
            design_flow_temp: 55.,
            target_flow_temp: 48.4375,
            temp_rm_prev: 22.425363294403255,
            variable_flow: true,
            temp_diff_emit_design: 10.,
            min_flow_rate: 0.05,
            max_flow_rate: 0.3,
        },
    ];

    #[rstest]
    fn test_buffer_loss(mut buffer_tank: BufferTank, simulation_time: SimulationTime) {
        // temperatures required to calculate buffer tank thermal losses
        let temp_ave_buffer = [
            42.13174056962889,
            39.15850036184071,
            47.90592627854362,
            46.396261107092855,
            44.393784168519964,
            41.113688636784985,
            38.05323632350841,
            34.94802876725388,
            37.658992822855524,
            34.3155583742916,
        ];
        let temp_rm_prev = [
            22.503414923272768,
            22.629952417028925,
            18.606533490587633,
            17.3014761483025,
            16.89603650204302,
            22.631060393299737,
            22.434635318905773,
            22.736995736582873,
            22.603763288379653,
            22.89241467168529,
        ];
        // Expected results
        let heat_loss_buffer_kwh = [
            0.01620674285529469,
            0.0063519154341823356,
            0.025287016057516827,
            0.010785181618173873,
            0.009663116173139813,
            0.0066316051216787795,
            0.01324052174653832,
            0.005063009401174876,
        ];
        let flow_temp_increase_due_to_buffer = [
            3.952751095382638,
            6.140725209053976,
            1.5784508035116716,
            3.8392108282420097,
            5.2146182138439485,
            7.521641387569073,
            7.370102324208936,
            6.569653721125626,
        ];
        let expected_thermal_losses = [
            0.015266475502721427,
            0.012855537290409166,
            0.022788416612854655,
            0.02262927719017028,
            0.02138713707392651,
            0.014375377522710748,
            0.012147800781357604,
            0.00949747013496634,
        ];
        let expected_internal_gains = [
            24.31011428294203,
            9.527873151273504,
            37.93052408627524,
            16.17777242726081,
            14.494674259709718,
            9.947407682518168,
            19.86078261980748,
            7.594514101762314,
        ];

        // Simulate updating buffer loss over time
        for (t_idx, _) in simulation_time.iter().enumerate() {
            let buffer_results = buffer_tank
                .calc_buffer_tank("test service", EMITTERS_DATA_FOR_BUFFER_TANK[t_idx])
                .unwrap()
                .to_owned();
            let internal_gains = buffer_tank.internal_gains();
            let thermal_losses =
                buffer_tank.thermal_losses(temp_ave_buffer[t_idx], temp_rm_prev[t_idx]);

            assert_relative_eq!(
                buffer_results[t_idx].heat_loss_buffer_kwh,
                heat_loss_buffer_kwh[t_idx],
                max_relative = 1e-6
            );
            assert_relative_eq!(
                buffer_results[t_idx].flow_temp_increase_due_to_buffer,
                flow_temp_increase_due_to_buffer[t_idx],
                max_relative = 1e-6
            );
            assert_relative_eq!(
                thermal_losses,
                expected_thermal_losses[t_idx],
                max_relative = 1e-6
            );
            assert_relative_eq!(
                internal_gains,
                expected_internal_gains[t_idx],
                max_relative = 1e-6
            );
        }
    }

    #[rstest]
    fn test_update_buffer_loss(mut buffer_tank: BufferTank) {
        buffer_tank.update_buffer_loss(5.);

        assert_eq!(buffer_tank.buffer_loss(), 5.);
    }

    #[rstest]
    fn test_calc_buffer_tank_over_max_flow(mut buffer_tank: BufferTank) {
        let mut data = EMITTERS_DATA_FOR_BUFFER_TANK[0];
        data.max_flow_rate = 0.1;

        let expected = &BufferTankServiceResult {
            _service_name: "new_service_buffer_tank".into(),
            _power_req_from_buffer_tank: 6.325422354229758,
            _temp_emitter_req: 43.32561228292832,
            _buffer_emitter_circ_flow_rate: 15.,
            flow_temp_increase_due_to_buffer: 9.109608401991274,
            pump_power_at_flow_rate: 0.04,
            heat_loss_buffer_kwh: 0.01620674285529469,
        };

        let result = &buffer_tank.calc_buffer_tank("new_service", data).unwrap()[0];

        assert_eq!(expected, result);
    }

    #[rstest]
    fn test_calc_buffer_tank_over_variable_flow(mut buffer_tank: BufferTank) {
        let mut data = EMITTERS_DATA_FOR_BUFFER_TANK[0];
        data.variable_flow = false;

        let expected = &BufferTankServiceResult {
            _service_name: "new_service_buffer_tank".into(),
            _power_req_from_buffer_tank: 6.325422354229758,
            _temp_emitter_req: 43.32561228292832,
            _buffer_emitter_circ_flow_rate: 15.,
            flow_temp_increase_due_to_buffer: 24.26646570859991,
            pump_power_at_flow_rate: 0.04,
            heat_loss_buffer_kwh: 0.01620674285529469,
        };

        let result = &buffer_tank.calc_buffer_tank("new_service", data).unwrap()[0];

        assert_eq!(expected, result);
    }

    #[rstest]
    fn test_calc_buffer_tank_no_power_req(mut buffer_tank: BufferTank) {
        let mut data = EMITTERS_DATA_FOR_BUFFER_TANK[0];

        data.power_req_from_buffer_tank = 0.;

        let expected = &BufferTankServiceResult {
            _service_name: "new_service_buffer_tank".into(),
            _power_req_from_buffer_tank: 0.,
            _temp_emitter_req: 43.32561228292832,
            _buffer_emitter_circ_flow_rate: 15.,
            flow_temp_increase_due_to_buffer: 0.,
            pump_power_at_flow_rate: 0.,
            heat_loss_buffer_kwh: -0.0019354000314273378,
        };

        let result = &buffer_tank.calc_buffer_tank("new_service", data).unwrap()[0];

        assert_eq!(expected, result);
    }

    #[rstest]
    fn test_calc_buffer_tank_over_emitter_flow(mut buffer_tank: BufferTank) {
        let mut data = EMITTERS_DATA_FOR_BUFFER_TANK[0];
        data.variable_flow = false;
        data.max_flow_rate = 5.;
        data.min_flow_rate = 10.;

        assert!(buffer_tank.calc_buffer_tank("new_service", data).is_err());
    }

    #[rstest]
    fn test_calc_buffer_tank_detailed_results(mut buffer_tank: BufferTank) {
        let data = EMITTERS_DATA_FOR_BUFFER_TANK[0];

        buffer_tank.calc_buffer_tank("new_service", data).unwrap();

        assert!(buffer_tank.detailed_results.is_some());

        let results = buffer_tank.detailed_results.unwrap();

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].len(), 1);

        let expected = &BufferTankServiceResult {
            _service_name: "new_service_buffer_tank".into(),
            _power_req_from_buffer_tank: 6.325422354229758,
            _temp_emitter_req: 43.32561228292832,
            _buffer_emitter_circ_flow_rate: 15.,
            flow_temp_increase_due_to_buffer: 3.952751095382638,
            pump_power_at_flow_rate: 0.04,
            heat_loss_buffer_kwh: 0.01620674285529469,
        };
        let result = &results[0][0];

        assert_eq!(expected, result);
    }

    #[rstest]
    fn test_interpolate_exhaust_air_heat_pump_test_data() {
        let mut data_eahp = vec![
            HeatPumpTestDatumInput {
                air_flow_rate: Some(100.0),
                test_letter: TestLetter::A,
                capacity: 5.0,
                cop: 2.0,
                design_flow_temp: 55.,
                temp_outlet: 55.,
                temp_source: 20.,
                temp_test: -7.,
                eahp_mixed_ext_air_ratio: Some(0.25),
            },
            HeatPumpTestDatumInput {
                air_flow_rate: Some(200.0),
                test_letter: TestLetter::A,
                capacity: 6.0,
                cop: 2.5,
                design_flow_temp: 55.,
                temp_outlet: 55.,
                temp_source: 20.,
                temp_test: -7.,
                eahp_mixed_ext_air_ratio: Some(0.75),
            },
            HeatPumpTestDatumInput {
                air_flow_rate: Some(100.0),
                test_letter: TestLetter::B,
                capacity: 5.5,
                cop: 2.4,
                design_flow_temp: 35.,
                temp_outlet: 34.,
                temp_source: 20.,
                temp_test: 2.,
                eahp_mixed_ext_air_ratio: Some(0.25),
            },
            HeatPumpTestDatumInput {
                air_flow_rate: Some(200.0),
                test_letter: TestLetter::B,
                capacity: 6.0,
                cop: 3.0,
                design_flow_temp: 35.,
                temp_outlet: 34.,
                temp_source: 20.,
                temp_test: 2.,
                eahp_mixed_ext_air_ratio: Some(0.75),
            },
        ];
        let data_eahp_interpolated = vec![
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::A,
                capacity: 5.4,
                cop: 2.2,
                design_flow_temp: 55.,
                temp_outlet: 55.,
                temp_source: 20.,
                temp_test: -7.,
                eahp_mixed_ext_air_ratio: Some(0.45),
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::B,
                capacity: 5.7,
                cop: 2.64,
                design_flow_temp: 35.,
                temp_outlet: 34.,
                temp_source: 20.,
                temp_test: 2.,
                eahp_mixed_ext_air_ratio: Some(0.45),
            },
        ];

        let (lowest_air_flow_in_test_data, data_eahp_func_result) =
            interpolate_exhaust_air_heat_pump_test_data(
                140.0,
                &data_eahp,
                HeatPumpSourceType::ExhaustAirMixed,
            )
            .unwrap();

        assert_eq!(
            lowest_air_flow_in_test_data, 100.0,
            "incorrect lowest air flow rate identified"
        );
        assert_eq!(
            data_eahp_func_result, data_eahp_interpolated,
            "incorrect interpolation of exhaust air heat pump test data"
        );
        data_eahp[1].temp_source = 21.;
        assert!(interpolate_exhaust_air_heat_pump_test_data(
            140.,
            &data_eahp,
            HeatPumpSourceType::ExhaustAirMixed
        )
        .is_err());
    }

    #[test]
    fn test_interpolate_exhaust_air_heat_pump_test_data_error() {
        assert!(interpolate_exhaust_air_heat_pump_test_data(
            140.,
            &vec![],
            HeatPumpSourceType::ExhaustAirMixed,
        )
        .is_err());
    }

    // NOTE from Python:
    // Before defining the code to run the tests, we define the data to be parsed
    // and the sorted/processed data structure it should be transformed into. Note
    // that the data for design flow temp of 55 has an extra record (test letter F2)
    // to test that the code can handle more than 2 records with the same temp_test
    // value properly. This probably won't occur in practice.
    #[fixture]
    pub fn data_unsorted() -> Vec<HeatPumpTestDatum> {
        vec![
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::A,
                capacity: 8.4,
                cop: 4.6,
                design_flow_temp: 35.,
                temp_outlet: 34.,
                temp_source: 0.,
                temp_test: -7.,
                eahp_mixed_ext_air_ratio: None,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::B,
                capacity: 8.3,
                cop: 4.9,
                design_flow_temp: 35.,
                temp_outlet: 30.,
                temp_source: 0.,
                temp_test: 2.,
                eahp_mixed_ext_air_ratio: None,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::C,
                capacity: 8.3,
                cop: 5.1,
                design_flow_temp: 35.,
                temp_outlet: 27.,
                temp_source: 0.,
                temp_test: 7.,
                eahp_mixed_ext_air_ratio: None,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::D,
                capacity: 8.2,
                cop: 5.4,
                design_flow_temp: 35.,
                temp_outlet: 24.,
                temp_source: 0.,
                temp_test: 12.,
                eahp_mixed_ext_air_ratio: None,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::F,
                capacity: 8.4,
                cop: 4.6,
                design_flow_temp: 35.,
                temp_outlet: 34.,
                temp_source: 0.,
                temp_test: -7.,
                eahp_mixed_ext_air_ratio: None,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::A,
                capacity: 8.8,
                cop: 3.2,
                design_flow_temp: 55.,
                temp_outlet: 52.,
                temp_source: 0.,
                temp_test: -7.,
                eahp_mixed_ext_air_ratio: None,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::B,
                capacity: 8.6,
                cop: 3.6,
                design_flow_temp: 55.,
                temp_outlet: 42.,
                temp_source: 0.,
                temp_test: 2.,
                eahp_mixed_ext_air_ratio: None,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::C,
                capacity: 8.5,
                cop: 3.9,
                design_flow_temp: 55.,
                temp_outlet: 36.,
                temp_source: 0.,
                temp_test: 7.,
                eahp_mixed_ext_air_ratio: None,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::D,
                capacity: 8.5,
                cop: 4.3,
                design_flow_temp: 55.,
                temp_outlet: 30.,
                temp_source: 0.,
                temp_test: 12.,
                eahp_mixed_ext_air_ratio: None,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::F,
                capacity: 8.8,
                cop: 3.2,
                design_flow_temp: 55.,
                temp_outlet: 52.,
                temp_source: 0.,
                temp_test: -7.,
                eahp_mixed_ext_air_ratio: None,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: TestLetter::F,
                capacity: 8.8,
                cop: 3.2,
                design_flow_temp: 55.,
                temp_outlet: 52.,
                temp_source: 0.,
                temp_test: -7.,
                eahp_mixed_ext_air_ratio: None,
            },
        ]
    }

    #[fixture]
    pub fn data_sorted() -> HashMap<OrderedFloat<f64>, Vec<CompleteHeatPumpTestDatum>> {
        let mut data: HashMap<_, _> = Default::default();
        data.insert(
            OrderedFloat(35.),
            vec![
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::A,
                    capacity: 8.4,
                    carnot_cop: 9.033823529411764,
                    cop: 4.6,
                    design_flow_temp: 35.,
                    exergetic_eff: 0.5091974605241738,
                    temp_outlet: 34.,
                    temp_source: 0.,
                    temp_test: -7.,
                    theoretical_load_ratio: 1.0,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::F,
                    capacity: 8.4,
                    carnot_cop: 9.033823529438331,
                    cop: 4.6,
                    design_flow_temp: 35.,
                    exergetic_eff: 0.5091974605226763,
                    temp_outlet: 34.,
                    temp_source: 0.0000000001,
                    temp_test: -6.9999999999,
                    theoretical_load_ratio: 1.0000000000040385,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::B,
                    capacity: 8.3,
                    carnot_cop: 10.104999999999999,
                    cop: 4.9,
                    design_flow_temp: 35.,
                    exergetic_eff: 0.48490846115784275,
                    temp_outlet: 30.,
                    temp_source: 0.,
                    temp_test: 2.,
                    theoretical_load_ratio: 1.1634388356892613,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::C,
                    capacity: 8.3,
                    carnot_cop: 11.116666666666665,
                    cop: 5.1,
                    design_flow_temp: 35.,
                    exergetic_eff: 0.4587706146926537,
                    temp_outlet: 27.,
                    temp_source: 0.,
                    temp_test: 7.,
                    theoretical_load_ratio: 1.3186802349509577,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::D,
                    capacity: 8.2,
                    carnot_cop: 12.38125,
                    cop: 5.4,
                    design_flow_temp: 35.,
                    exergetic_eff: 0.43614336193841496,
                    temp_outlet: 24.,
                    temp_source: 0.,
                    temp_test: 12.,
                    theoretical_load_ratio: 1.513621351820552,
                },
            ],
        );
        data.insert(
            OrderedFloat(55.),
            vec![
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::A,
                    capacity: 8.8,
                    carnot_cop: 6.252884615384615,
                    cop: 3.2,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.5117638013224666,
                    temp_outlet: 52.,
                    temp_source: 0.,
                    temp_test: -7.,
                    theoretical_load_ratio: 1.0,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::F,
                    capacity: 8.8,
                    carnot_cop: 6.252884615396638,
                    cop: 3.2,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.5117638013214826,
                    temp_outlet: 52.,
                    temp_source: 0.0000000001,
                    temp_test: -6.9999999999,
                    theoretical_load_ratio: 1.0000000000030207,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::F,
                    capacity: 8.8,
                    carnot_cop: 6.252884615408662,
                    cop: 3.2,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.5117638013204985,
                    temp_outlet: 52.,
                    temp_source: 0.0000000002,
                    temp_test: -6.9999999998,
                    theoretical_load_ratio: 1.0000000000060418,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::B,
                    capacity: 8.6,
                    carnot_cop: 7.503571428571428,
                    cop: 3.6,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.4797715373631604,
                    temp_outlet: 42.,
                    temp_source: 0.,
                    temp_test: 2.,
                    theoretical_load_ratio: 1.3179136223360988,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::C,
                    capacity: 8.5,
                    carnot_cop: 8.587499999999999,
                    cop: 3.9,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.4541484716157206,
                    temp_outlet: 36.,
                    temp_source: 0.,
                    temp_test: 7.,
                    theoretical_load_ratio: 1.5978273764295179,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: TestLetter::D,
                    capacity: 8.5,
                    carnot_cop: 10.104999999999999,
                    cop: 4.3,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.4255319148936171,
                    temp_outlet: 30.,
                    temp_source: 0.,
                    temp_test: 12.,
                    theoretical_load_ratio: 1.9940427298329144,
                },
            ],
        );

        data
    }

    #[fixture]
    pub fn test_data(data_unsorted: Vec<HeatPumpTestDatum>) -> HeatPumpTestData {
        HeatPumpTestData::new(data_unsorted).unwrap()
    }

    #[rstest]
    // In Python this test is called `test_init`
    pub fn should_have_constructed_internal_data_structures(
        test_data: HeatPumpTestData,
        data_sorted: HashMap<OrderedFloat<f64>, Vec<CompleteHeatPumpTestDatum>>,
    ) {
        assert_eq!(
            test_data.design_flow_temps,
            vec![OrderedFloat(35.), OrderedFloat(55.)],
            "list of design flow temps populated incorrectly"
        );
        assert_eq!(
            test_data
                .test_data
                .iter()
                .map(|(key, val)| {
                    (
                        key,
                        val.iter()
                            .map(|datum| datum.round_by_precision(1e7))
                            .collect(),
                    )
                })
                .collect::<HashMap<&OrderedFloat<f64>, Vec<CompleteHeatPumpTestDatum>>>(),
            data_sorted
                .iter()
                .map(|(key, val)| {
                    (
                        key,
                        val.iter()
                            .map(|datum| datum.round_by_precision(1e7))
                            .collect(),
                    )
                })
                .collect::<HashMap<&OrderedFloat<f64>, Vec<CompleteHeatPumpTestDatum>>>(),
            "list of test data records populated incorrectly"
        );
    }

    #[test]
    fn test_init_no_temps() {
        assert!(HeatPumpTestData::new(vec![]).is_err());
    }

    #[rstest]
    fn test_init_expected_five_letters(data_unsorted: Vec<HeatPumpTestDatum>) {
        // Test that init throws if there aren't 5 records for each design flow temperature
        let mut data = data_unsorted.clone();
        data.truncate(9);

        assert!(HeatPumpTestData::new(data).is_err());
    }

    #[rstest]
    fn test_init_maximum_design_flow() {
        // Test that initialising with more than 4 design flow temperatures errors
        let design_flow_temps = [35., 45., 55., 20., 65.];

        let base_values = [
            (TestLetter::A, 8.4, 4.6, 34., 0., -7.),
            (TestLetter::B, 8.3, 4.9, 30., 0., 2.),
            (TestLetter::C, 8.3, 5.1, 27., 0., 7.),
            (TestLetter::D, 8.2, 5.4, 24., 0., 12.),
            (TestLetter::F, 8.4, 4.6, 34., 0., -7.),
        ];

        let data: Vec<HeatPumpTestDatum> = design_flow_temps
            .iter()
            .flat_map(|design_flow_temp| {
                base_values.iter().map(
                    |&(test_letter, capacity, cop, temp_outlet, temp_source, temp_test)| {
                        HeatPumpTestDatum {
                            air_flow_rate: None,
                            test_letter,
                            capacity,
                            cop,
                            design_flow_temp: *design_flow_temp,
                            temp_outlet,
                            temp_source,
                            temp_test,
                            eahp_mixed_ext_air_ratio: None,
                        }
                    },
                )
            })
            .collect();

        assert!(HeatPumpTestData::new(data).is_err());
    }

    #[rstest]
    fn test_init_expected_letter(data_unsorted: Vec<HeatPumpTestDatum>) {
        // Test that initialisation errors when test letter F is missing
        let mut data: Vec<HeatPumpTestDatum> = data_unsorted.clone();
        data.truncate(9);
        data.push(HeatPumpTestDatum {
            air_flow_rate: None,
            test_letter: TestLetter::A, // In Python this is "E"
            capacity: 8.5,
            cop: 4.3,
            design_flow_temp: 55.,
            temp_outlet: 30.,
            temp_source: 0.,
            temp_test: 12.,
            eahp_mixed_ext_air_ratio: None,
        });

        assert!(HeatPumpTestData::new(data).is_err());
    }

    #[rstest]
    fn test_init_four_distinct_records(data_unsorted: Vec<HeatPumpTestDatum>) {
        // Test that initialisation errors if there aren't four distinct records for each design flow temperature
        let mut data = data_unsorted.clone();
        data.truncate(5);

        data[0].temp_test = 1.;
        data[1].temp_test = 1.;
        data[2].temp_test = 1.;
        data[3].temp_test = 1.;
        data[4].temp_test = 1.;

        assert!(HeatPumpTestData::new(data).is_err());
    }

    #[fixture]
    fn expected_init_regression_coeffs() -> HashMap<OrderedFloat<f64>, Vec<f64>> {
        let mut expected: HashMap<_, _> = Default::default();
        expected.insert(
            OrderedFloat(35.),
            vec![
                4.810017281274474,
                0.03677543129969712,
                0.0009914765238219557,
            ],
        );
        expected.insert(
            OrderedFloat(55.),
            vec![
                3.4857982546529747,
                0.050636568790103545,
                0.0014104955583514216,
            ],
        );
        expected
    }

    #[rstest]
    fn test_init_regression_coeffs(
        test_data: HeatPumpTestData,
        expected_init_regression_coeffs: HashMap<OrderedFloat<f64>, Vec<f64>>,
    ) {
        assert_eq!(
            test_data
                .regression_coeffs
                .iter()
                .map(|(key, val)| (key, round_each_by_precision(val, 1e7)))
                .collect::<HashMap<_, _>>(),
            expected_init_regression_coeffs
                .iter()
                .map(|(key, val)| (key, round_each_by_precision(val, 1e7)))
                .collect::<HashMap<_, _>>(),
            "list of regression coefficients populated incorrectly"
        );
    }

    #[rstest]
    fn test_average_capacity(test_data: HeatPumpTestData) {
        let results = [8.3, 8.375, 8.45, 8.525, 8.6];
        for (i, flow_temp) in [35., 40., 45., 50., 55.].iter().enumerate() {
            assert_ulps_eq!(
                test_data
                    .average_capacity(celsius_to_kelvin(*flow_temp).unwrap())
                    .unwrap(),
                results[i],
            )
        }
    }

    #[rstest]
    fn test_temp_spread_test_conditions(test_data: HeatPumpTestData) {
        let results = [5.0, 5.75, 6.5, 7.25, 8.0];
        for (i, flow_temp) in [35., 40., 45., 50., 55.].iter().enumerate() {
            assert_eq!(
                test_data
                    .temp_spread_test_conditions(celsius_to_kelvin(*flow_temp).unwrap())
                    .unwrap(),
                results[i],
                "incorrect temp spread at test conditions returned"
            );
        }
    }

    #[rstest]
    fn test_find_test_record_index(test_data: HeatPumpTestData) {
        assert!(test_data.find_test_record_index("x", 55.).is_err());
    }

    #[fixture]
    pub fn carnot_cop_cases() -> Vec<(f64, &'static str, f64)> {
        vec![
            (35., "cld", 9.033823529411764),
            (40., "cld", 8.338588800904978),
            (45., "cld", 7.643354072398189),
            (50., "cld", 6.948119343891403),
            (55., "cld", 6.252884615384615),
            (45., "A", 7.643354072398189),
            (45., "B", 8.804285714285713),
            (45., "C", 9.852083333333333),
            (45., "D", 11.243125),
            (45., "F", 7.643354072417485),
        ]
        // TODO (from Python) Note that the result above for condition F is different to
        // that for condition A, despite the source and outlet temps in
        // the inputs being the same for both, because of the adjustment
        // to the source temp applied in the HeatPumpTestData __init__
        // function when duplicate records are found. This may not be
        // the desired behaviour (see the TODO comment in that function)
        // but in that case the problem is not with the function that
        // is being tested here, so for now we set the result so that
        // the test passes.
    }

    #[rstest]
    pub fn test_carnot_cop_at_test_condition(
        test_data: HeatPumpTestData,
        carnot_cop_cases: Vec<(f64, &str, f64)>,
    ) {
        // TODO (from Python) Test conditions other than just coldest
        for (flow_temp, test_condition, result) in carnot_cop_cases {
            assert_ulps_eq!(
                test_data
                    .carnot_cop_at_test_condition(
                        test_condition,
                        celsius_to_kelvin(flow_temp).unwrap(),
                    )
                    .unwrap(),
                result,
            );
        }
    }

    #[fixture]
    pub fn outlet_temp_cases() -> Vec<(f64, &'static str, f64)> {
        // TODO (from Python) Test conditions other than just coldest
        vec![
            (35., "cld", 307.15),
            (40., "cld", 311.65),
            (45., "cld", 316.15),
            (50., "cld", 320.65),
            (55., "cld", 325.15),
            (45., "A", 316.15),
            (45., "B", 309.15),
            (45., "C", 304.65),
            (45., "D", 300.15),
            (45., "F", 316.15),
        ]
    }

    #[rstest]
    pub fn test_outlet_temp_at_test_condition(
        test_data: HeatPumpTestData,
        outlet_temp_cases: Vec<(f64, &'static str, f64)>,
    ) {
        for (flow_temp, test_condition, result) in outlet_temp_cases {
            assert_ulps_eq!(
                test_data
                    .outlet_temp_at_test_condition(
                        test_condition,
                        celsius_to_kelvin(flow_temp).unwrap(),
                    )
                    .unwrap(),
                result,
            );
        }
    }

    #[fixture]
    pub fn source_temp_cases() -> Vec<(f64, &'static str, f64)> {
        // TODO (from Python) Test conditions other than just coldest
        vec![
            (35., "cld", 273.15),
            (40., "cld", 273.15),
            (45., "cld", 273.15),
            (50., "cld", 273.15),
            (55., "cld", 273.15),
            (45., "A", 273.15),
            (45., "B", 273.15),
            (45., "C", 273.15),
            (45., "D", 273.15),
            (45., "F", 273.15000000009996),
        ]

        // TODO (from Python) Note that the result above for condition F is different to
        // that for condition A, despite the source and outlet temps in
        // the inputs being the same for both, because of the adjustment
        // to the source temp applied in the HeatPumpTestData __init__
        // function when duplicate records are found. This may not be
        // the desired behaviour (see the TODO comment in that function)
        // but in that case the problem is not with the function that
        // is being tested here, so for now we set the result so that
        // the test passes.
    }

    #[rstest]
    pub fn test_source_temp_at_test_condition(
        test_data: HeatPumpTestData,
        source_temp_cases: Vec<(f64, &'static str, f64)>,
    ) {
        for (flow_temp, test_condition, result) in source_temp_cases {
            assert_ulps_eq!(
                test_data
                    .source_temp_at_test_condition(
                        test_condition,
                        celsius_to_kelvin(flow_temp).unwrap(),
                    )
                    .unwrap(),
                result,
            );
        }
    }

    #[fixture]
    pub fn capacity_cases() -> Vec<(f64, &'static str, f64)> {
        vec![
            (35., "cld", 8.4),
            (40., "cld", 8.5),
            (45., "cld", 8.6),
            (50., "cld", 8.7),
            (55., "cld", 8.8),
            (45., "A", 8.6),
            (45., "B", 8.45),
            (45., "C", 8.4),
            (45., "D", 8.35),
            (45., "F", 8.6),
        ]
    }

    #[rstest]
    pub fn test_capacity_at_test_condition(
        test_data: HeatPumpTestData,
        capacity_cases: Vec<(f64, &'static str, f64)>,
    ) {
        for (flow_temp, test_condition, result) in capacity_cases {
            assert_ulps_eq!(
                test_data
                    .capacity_at_test_condition(
                        test_condition,
                        celsius_to_kelvin(flow_temp).unwrap()
                    )
                    .unwrap(),
                result,
            );
        }
    }

    #[fixture]
    pub fn lr_op_cond_cases() -> Vec<[f64; 4]> {
        vec![
            [35.0, 283.15, 12.326, 2.579601322115558],
            [40.0, 293.15, 15.6575, 3.464991021429274],
            [45.0, 278.15, 7.95375, 1.4337736417538745],
            [50.0, 288.15, 9.23285714285714, 1.765821306168971],
            [55.0, 273.15, 5.96636363636364, 0.9282465711109492],
        ]
    }

    fn design_flow_temp_op_cond_k(design_flow_temp_op_cond: f64) -> f64 {
        design_flow_temp_op_cond + 273.15
    }

    #[rstest]
    // In Python this test is called `test_lr_op_cond`
    pub fn test_load_ratio_under_operational_conditions(
        test_data: HeatPumpTestData,
        lr_op_cond_cases: Vec<[f64; 4]>,
    ) {
        for [flow_temp, temp_source, carnot_cop_op_cond, result] in lr_op_cond_cases {
            assert_relative_eq!(
                test_data
                    .load_ratio_at_operating_conditions(
                        celsius_to_kelvin(flow_temp).unwrap(),
                        temp_source,
                        carnot_cop_op_cond,
                        design_flow_temp_op_cond_k(45.),
                    )
                    .unwrap(),
                result,
                max_relative = 1e-7
            );
        }
    }

    #[rstest]
    pub fn test_lr_eff_either_side_of_op_cond(test_data: HeatPumpTestData) {
        let results_lr_below = [
            1.1634388356892613,
            1.0000000000060418,
            1.3186802349509577,
            1.3179136223360988,
            1.3186802349509577,
            1.5978273764295179,
        ];
        let results_lr_above = [
            1.3186802349509577,
            1.3179136223360988,
            1.513621351820552,
            1.5978273764295179,
            1.513621351820552,
            1.9940427298329144,
        ];
        let results_eff_below = [
            0.48490846115784275,
            0.5117638013204985,
            0.4587706146926537,
            0.4797715373631604,
            0.4587706146926537,
            0.4541484716157206,
        ];
        let results_eff_above = [
            0.4587706146926537,
            0.4797715373631604,
            0.43614336193841496,
            0.4541484716157206,
            0.43614336193841496,
            0.4255319148936171,
        ];

        let mut i = 0;
        for exergy_lr_op_cond in [1.2, 1.4, 1.6] {
            for flow_temp in [35., 55.] {
                let (lr_below, lr_above, eff_below, eff_above) = test_data
                    .lr_eff_either_side_of_op_cond(&OrderedFloat(flow_temp), exergy_lr_op_cond)
                    .unwrap();
                assert_relative_eq!(lr_below, results_lr_below[i], max_relative = 1e-7);
                assert_ulps_eq!(lr_above, results_lr_above[i],);
                assert_relative_eq!(eff_below, results_eff_below[i], max_relative = 1e-7);
                assert_ulps_eq!(eff_above, results_eff_above[i],);
                i += 1;
            }
        }

        // first load ratio in the test data cannot be greater than the load ratio at operating conditions
        assert_eq!(
            test_data
                .lr_eff_either_side_of_op_cond(&OrderedFloat(55.), 0.9)
                .unwrap(),
            (
                1.,
                1.0000000000030207,
                0.5117638013224666,
                0.5117638013214826
            )
        );
    }

    #[rstest]
    pub fn test_cop_op_cond_if_not_air_source(test_data: HeatPumpTestData) {
        let results = [
            6.717150897917277,
            8.121968575142535,
            4.60977003063163,
            5.748091483057522,
            3.579953030724001,
        ];
        let temp_cases = [
            [8.0, 0.00, 283.15, 308.15],
            [7.0, -5.0, 293.15, 313.15],
            [6.0, 5.00, 278.15, 318.15],
            [5.0, 10.0, 288.15, 323.15],
            [4.0, 7.50, 273.15, 328.15],
        ];

        for (i, [temp_diff_limit_low, temp_ext, temp_source, temp_output]) in
            temp_cases.iter().enumerate()
        {
            assert_relative_eq!(
                test_data
                    .cop_op_cond_if_not_air_source(
                        *temp_diff_limit_low,
                        *temp_ext,
                        *temp_source,
                        *temp_output,
                        design_flow_temp_op_cond_k(45.),
                    )
                    .unwrap(),
                results[i],
                max_relative = 1e-7
            );
        }
    }

    #[rstest]
    pub fn test_capacity_op_cond_var_flow_or_source_temp(test_data: HeatPumpTestData) {
        let results = [
            10.390871958470573,
            8.241818181818182,
            8.95014809894768,
            9.495814072084585,
            8.830454545454547,
        ];
        let test_cases = [
            (true, 283.15, 308.15),
            (false, 293.15, 313.15),
            (true, 278.15, 318.15),
            (true, 288.15, 323.15),
            (false, 273.15, 328.15),
        ];
        for (i, (mod_ctrl, temp_source, temp_output)) in test_cases.iter().enumerate() {
            assert_relative_eq!(
                test_data
                    .capacity_op_cond_var_flow_or_source_temp(
                        *temp_output,
                        *temp_source,
                        *mod_ctrl,
                        design_flow_temp_op_cond_k(45.)
                    )
                    .unwrap(),
                results[i],
                max_relative = 1e-7
            );
        }
    }

    #[rstest]
    pub fn test_temp_spread_correction(test_data: HeatPumpTestData) {
        let results = [
            1.0872913992297817,
            1.0698529411764706,
            1.05822498586772,
            1.0499171499585749,
            1.043684710351377,
        ];
        let temp_source = 275.15;
        let temp_diff_evaporator = -15.0;
        let temp_diff_condenser = 5.0;
        let temp_spread_emitter = 10.0;

        for (i, temp_output) in [308.15, 313.15, 318.15, 323.15, 328.15].iter().enumerate() {
            assert_relative_eq!(
                test_data
                    .temp_spread_correction(
                        temp_source,
                        *temp_output,
                        temp_diff_evaporator,
                        temp_diff_condenser,
                        temp_spread_emitter,
                        design_flow_temp_op_cond_k(45.)
                    )
                    .unwrap(),
                results[i],
                max_relative = 1e-7
            );
        }
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }

    fn round_each_by_precision(src: &[f64], precision: f64) -> Vec<f64> {
        src.iter()
            .map(|num| round_by_precision(*num, precision))
            .collect::<Vec<f64>>()
    }

    impl CompleteHeatPumpTestDatum {
        pub fn round_by_precision(&self, precision: f64) -> CompleteHeatPumpTestDatum {
            CompleteHeatPumpTestDatum {
                air_flow_rate: self.air_flow_rate,
                test_letter: self.test_letter,
                capacity: self.capacity,
                cop: self.cop,
                design_flow_temp: self.design_flow_temp,
                temp_outlet: self.temp_outlet,
                temp_source: round_by_precision(self.temp_source, precision),
                temp_test: round_by_precision(self.temp_test, precision),
                carnot_cop: round_by_precision(self.carnot_cop, precision),
                exergetic_eff: round_by_precision(self.exergetic_eff, precision),
                theoretical_load_ratio: round_by_precision(self.theoretical_load_ratio, precision),
            }
        }
    }

    #[rstest]
    // In Python this test is called `test_from_string` (inside the `TestSourceType` class)
    pub fn test_source_type_from_string() {
        let types = [
            ("Ground", HeatPumpSourceType::Ground),
            ("OutsideAir", HeatPumpSourceType::OutsideAir),
            ("ExhaustAirMEV", HeatPumpSourceType::ExhaustAirMEV),
            ("ExhaustAirMVHR", HeatPumpSourceType::ExhaustAirMVHR),
            ("ExhaustAirMixed", HeatPumpSourceType::ExhaustAirMixed),
            ("WaterGround", HeatPumpSourceType::WaterGround),
            ("WaterSurface", HeatPumpSourceType::WaterSurface),
            ("HeatNetwork", HeatPumpSourceType::HeatNetwork),
        ];
        for (strval, source_type) in types {
            assert_eq!(
                serde_json::from_str::<HeatPumpSourceType>(format!("\"{strval}\"").as_str())
                    .unwrap(),
                source_type
            );
        }
    }

    #[rstest]
    // In Python this test is called `test_is_exhaust_air` (inside the `TestSourceType` class)
    pub fn test_source_type_is_exhaust_air() {
        let exhaust_air_cases = [
            (HeatPumpSourceType::Ground, false),
            (HeatPumpSourceType::OutsideAir, false),
            (HeatPumpSourceType::ExhaustAirMEV, true),
            (HeatPumpSourceType::ExhaustAirMVHR, true),
            (HeatPumpSourceType::ExhaustAirMixed, true),
            (HeatPumpSourceType::WaterGround, false),
            (HeatPumpSourceType::WaterSurface, false),
            (HeatPumpSourceType::HeatNetwork, false),
        ];
        for (source_type, result) in exhaust_air_cases {
            assert_eq!(
                source_type.is_exhaust_air(),
                result,
                "incorrectly identified whether or not source type is exhaust air"
            );
        }
    }

    #[rstest]
    pub fn test_source_fluid_is_air() {
        let fluid_is_air_cases = [
            (HeatPumpSourceType::Ground, false),
            (HeatPumpSourceType::OutsideAir, true),
            (HeatPumpSourceType::ExhaustAirMEV, true),
            (HeatPumpSourceType::ExhaustAirMVHR, true),
            (HeatPumpSourceType::ExhaustAirMixed, true),
            (HeatPumpSourceType::WaterGround, false),
            (HeatPumpSourceType::WaterSurface, false),
            (HeatPumpSourceType::HeatNetwork, false),
        ];
        for (source_type, result) in fluid_is_air_cases {
            assert_eq!(
                source_type.source_fluid_is_air(),
                result,
                "incorrectly identified whether or not source fluid is air"
            );
        }
    }

    #[rstest]
    pub fn test_source_fluid_is_water() {
        let fluid_is_water_cases = [
            (HeatPumpSourceType::Ground, true),
            (HeatPumpSourceType::OutsideAir, false),
            (HeatPumpSourceType::ExhaustAirMEV, false),
            (HeatPumpSourceType::ExhaustAirMVHR, false),
            (HeatPumpSourceType::ExhaustAirMixed, false),
            (HeatPumpSourceType::WaterGround, true),
            (HeatPumpSourceType::WaterSurface, true),
            (HeatPumpSourceType::HeatNetwork, true),
        ];
        for (source_type, result) in fluid_is_water_cases {
            assert_eq!(
                source_type.source_fluid_is_water(),
                result,
                "incorrectly identified whether or not source fluid is water"
            );
        }
    }

    // TestHeatPumpService
    #[rstest]
    fn test_is_on_for_service_water(
        heat_pump_service_water: HeatPumpServiceWater,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        assert!(
            heat_pump_service_water.is_on(simulation_time_for_heat_pump.iter().current_iteration())
        )
    }

    #[rstest]
    fn test_is_on_for_service_space(
        heat_pump_service_space: HeatPumpServiceSpace,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        assert!(
            heat_pump_service_space.is_on(simulation_time_for_heat_pump.iter().current_iteration())
        )
    }

    #[rstest]
    fn test_is_on_for_service_space_warm_air(
        mut heat_pump_service_space_warm_air: HeatPumpServiceSpaceWarmAir,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        heat_pump_service_space_warm_air.control =
            Arc::new(create_setpoint_time_control(vec![Some(10.)]));

        assert!(heat_pump_service_space_warm_air
            .is_on(simulation_time_for_heat_pump.iter().current_iteration()))
    }

    // TestHeatPumpServiceWater
    #[fixture]
    fn heat_pump_service_water(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) -> HeatPumpServiceWater {
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
        let control_min = create_setpoint_time_control(vec![Some(10.)]);
        let control_max = create_setpoint_time_control(vec![Some(20.)]);
        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(
            ColdWaterSource::new(vec![1.0, 1.2], 0, simulation_time_for_heat_pump.step).into(),
        ); // Python uses a mock for cold_feed

        HeatPumpServiceWater::new(
            Arc::new(Mutex::new(heat_pump)),
            "new_service".into(),
            30.,
            Arc::new(cold_feed),
            Arc::new(control_min),
            Arc::new(control_max),
            None,
        )
    }

    #[rstest]
    fn test_setpnt_for_service_water(
        heat_pump_service_water: HeatPumpServiceWater,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let (setpnt_min, setpnt_max) = heat_pump_service_water
            .setpnt(simulation_time_for_heat_pump.iter().current_iteration())
            .unwrap();

        assert_eq!(setpnt_min, Some(10.));
        assert_eq!(setpnt_max, Some(20.));
    }

    #[rstest]
    fn test_setpnt_errors_for_service_water(
        mut heat_pump_service_water: HeatPumpServiceWater,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        heat_pump_service_water.control_min = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(true)],
            0,
            1.,
        )));

        assert!(heat_pump_service_water
            .setpnt(simulation_time_for_heat_pump.iter().current_iteration())
            .is_err());

        heat_pump_service_water.control_min =
            Arc::new(create_setpoint_time_control(vec![Some(10.)]));

        heat_pump_service_water.control_max = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(true)],
            0,
            1.,
        )));

        assert!(heat_pump_service_water
            .setpnt(simulation_time_for_heat_pump.iter().current_iteration())
            .is_err());
    }

    // skipping Python's test_energy_output_max due to mocking
    // skipping Python's test_energy_output_max_no_temp_return due to mocking
    // skipping Python's test_demand_energy due to mocking
    // skipping Python's test_demand_energy_no_temp_return due to mocking
    // skipping Python's test_demand_energy_no_temp_flow due to mocking

    // TestHeatPumpServiceSpace
    #[fixture]
    fn heat_pump_service_space(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) -> HeatPumpServiceSpace {
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
        let control = create_setpoint_time_control(vec![Some(20.), None]);

        HeatPumpServiceSpace::new(
            Arc::new(Mutex::new(heat_pump)),
            "new_service".into(),
            HeatPumpEmitterType::RadiatorsUfh,
            50.,
            Default::default(),
            55.,
            Arc::new(control),
            100.,
            None,
        )
    }

    #[rstest]
    fn test_setpnt_for_service_space(
        heat_pump_service_space: HeatPumpServiceSpace,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        assert_eq!(
            heat_pump_service_space
                .temp_setpnt(&simulation_time_for_heat_pump.iter().current_iteration())
                .unwrap(),
            Some(20.)
        );
    }

    #[rstest]
    fn test_in_required_period(
        heat_pump_service_space: HeatPumpServiceSpace,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            assert_eq!(
                heat_pump_service_space.in_required_period(&t_it).unwrap(),
                [Some(true), Some(false)][t_idx]
            );
        }
    }

    #[rstest]
    fn test_control_errors(
        mut heat_pump_service_space: HeatPumpServiceSpace,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        heat_pump_service_space.control = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(true)],
            0,
            1.,
        )));

        assert!(heat_pump_service_space
            .temp_setpnt(&simulation_time_for_heat_pump.iter().current_iteration())
            .is_err());
        assert!(heat_pump_service_space
            .in_required_period(&simulation_time_for_heat_pump.iter().current_iteration())
            .is_err());
    }

    // skipping Python's test_energy_output_max due to mocking

    #[rstest]
    fn test_energy_output_max_off_for_service_space(
        mut heat_pump_service_space: HeatPumpServiceSpace,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let control = create_setpoint_time_control(vec![None]);
        heat_pump_service_space.control = Arc::new(control);

        assert_eq!(
            heat_pump_service_space
                .energy_output_max(
                    50.,
                    40.,
                    None,
                    None,
                    simulation_time_for_heat_pump.iter().current_iteration()
                )
                .unwrap()
                .0,
            0.
        );
    }

    // skipping Python's test_demand_energy due to mocking
    // skipping Python's test_demand_energy_off due to mocking
    // skipping Python's test_running_time_throughput_factor due to mocking
    // skipping Python's test_running_time_throughput_factor_off due to mocking
    // skipping Python's test_temp_spread_correction_air due to mocking
    // skipping Python's test_temp_spread_correction_water due to mocking

    // skipping Python's test_temp_spread_correction_invalid as it's not possible to have an invalid source type in the Rust

    // TestHeatPumpServiceSpaceWarmAir
    #[fixture]
    fn heat_pump_service_space_warm_air(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) -> HeatPumpServiceSpaceWarmAir {
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
        let control = create_setpoint_time_control(vec![]);

        HeatPumpServiceSpaceWarmAir::new(
            Arc::new(Mutex::new(heat_pump)),
            "new_service",
            10.,
            25.,
            Arc::new(control),
            10.,
            0.7,
            12.,
        )
    }

    #[rstest]
    fn test_energy_output_min(heat_pump_service_space_warm_air: HeatPumpServiceSpaceWarmAir) {
        assert_eq!(heat_pump_service_space_warm_air.energy_output_min(), 0.);
    }

    // skipping Python's test_demand_energy due to mocking
    // skipping Python's test_running_time_throughput_factor due to mocking

    #[rstest]
    fn test_frac_convective(heat_pump_service_space_warm_air: HeatPumpServiceSpaceWarmAir) {
        assert_eq!(heat_pump_service_space_warm_air.frac_convective(), 0.7);
    }

    #[rstest]
    // In Python this test is called `test_from_string` (inside the `TestSinkType` class)
    pub fn test_sink_type_from_string() {
        let sink_types = [
            ("Air", HeatPumpSinkType::Air),
            ("Water", HeatPumpSinkType::Water),
        ];
        for (strval, result) in sink_types {
            assert_eq!(
                serde_json::from_str::<HeatPumpSinkType>(format!("\"{strval}\"").as_str()).unwrap(),
                result,
                "incorrect source type mapped"
            );
        }
    }

    // The below tests correspond to the TestHeatPump class in Python

    #[fixture]
    fn simulation_time_for_heat_pump() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    fn external_conditions(simulation_time_for_heat_pump: SimulationTime) -> ExternalConditions {
        let simulation_time_iterator = simulation_time_for_heat_pump.iter();
        let wind_speeds = vec![3.7, 3.8];
        let wind_directions = vec![200., 220.];
        let air_temps = vec![0.0, 2.5];
        let diffuse_horizontal_radiations = vec![333., 610.];
        let direct_beam_radiations = vec![420., 750.];
        let shading_segments = vec![
            ShadingSegment {
                start: 180.,
                end: 135.,
                ..Default::default()
            },
            ShadingSegment {
                start: 135.,
                end: 90.,
                ..Default::default()
            },
        ]
        .into();
        ExternalConditions::new(
            &simulation_time_iterator,
            air_temps,
            wind_speeds,
            wind_directions,
            diffuse_horizontal_radiations,
            direct_beam_radiations,
            vec![0.2; 8760],
            51.42,
            -0.75,
            0,
            0,
            None,
            1.0,
            Some(1),
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            shading_segments,
        )
    }

    fn create_temp_internal_air_fn(canned_value: f64) -> TempInternalAirFn {
        Arc::new(move || canned_value)
    }

    #[fixture]
    fn energy_supply(simulation_time_for_heat_pump: SimulationTime) -> EnergySupply {
        EnergySupplyBuilder::new(
            FuelType::MainsGas,
            simulation_time_for_heat_pump.iter().total_steps(),
        )
        .build()
    }

    fn create_boiler(
        external_conditions: ExternalConditions,
        energy_supply: EnergySupply,
        simulation_time_for_heat_pump: SimulationTime,
        energy_supply_conn_name_auxiliary: &str,
    ) -> Boiler {
        let boiler_details = HeatSourceWetDetails::Boiler {
            energy_supply: "mains gas".into(),
            energy_supply_aux: "mains elec".into(),
            rated_power: 24.,
            efficiency_full_load: 0.891,
            efficiency_part_load: 0.991,
            boiler_location: HeatSourceLocation::Internal,
            modulation_load: 0.3,
            electricity_circ_pump: 0.06,
            electricity_part_load: 0.0131,
            electricity_full_load: 0.0388,
            electricity_standby: 0.0244,
        };

        let energy_supply: Arc<RwLock<EnergySupply>> = Arc::from(RwLock::from(energy_supply));

        let energy_supply_conn_aux = EnergySupplyConnection::new(
            energy_supply.clone(),
            energy_supply_conn_name_auxiliary.into(),
        );

        Boiler::new(
            boiler_details,
            energy_supply,
            energy_supply_conn_aux,
            Arc::new(external_conditions),
            simulation_time_for_heat_pump.step,
        )
        .unwrap()
    }

    fn create_heat_pump_input_from_json(
        backup_ctrl_type: Option<&str>,
        source_type: Option<&str>,
    ) -> HeatSourceWetDetails {
        let backup_ctrl_type = backup_ctrl_type.unwrap_or("TopUp");
        let source_type = source_type.unwrap_or("OutsideAir");
        let input = json!({
            "type": "HeatPump",
            "EnergySupply": "mains gas",
            "source_type": source_type,
            "sink_type": "Water",
            "backup_ctrl_type": backup_ctrl_type,
            "time_delay_backup": 1,
            "modulating_control": true,
            "min_modulation_rate_35": 0.35,
            "min_modulation_rate_55": 0.4,
            "time_constant_onoff_operation": 140,
            "temp_return_feed_max": 70,
            "temp_lower_operating_limit": -5,
            "min_temp_diff_flow_return_for_hp_to_operate": 0,
            "var_flow_temp_ctrl_during_test": false,
            "power_heating_circ_pump": 0.015,
            "power_source_circ_pump": 0.01,
            "power_standby": 0.015,
            "power_crankcase_heater": 0.01,
            "power_off": 0.015,
            "power_max_backup": 3,
            "test_data_EN14825": [
                {
                    "test_letter": "A",
                    "capacity": 8.4,
                    "cop": 4.6,
                    "design_flow_temp": 35,
                    "temp_outlet": 34,
                    "temp_source": 0,
                    "temp_test": -7
                },
                {
                    "test_letter": "B",
                    "capacity": 8.3,
                    "cop": 4.9,
                    "design_flow_temp": 35,
                    "temp_outlet": 30,
                    "temp_source": 0,
                    "temp_test": 2
                },
                {
                    "test_letter": "C",
                    "capacity": 8.3,
                    "cop": 5.1,
                    "design_flow_temp": 35,
                    "temp_outlet": 27,
                    "temp_source": 0,
                    "temp_test": 7
                },
                {
                    "test_letter": "D",
                    "capacity": 8.2,
                    "cop": 5.4,
                    "design_flow_temp": 35,
                    "temp_outlet": 24,
                    "temp_source": 0,
                    "temp_test": 12
                },
                {
                    "test_letter": "F",
                    "capacity": 8.4,
                    "cop": 4.6,
                    "design_flow_temp": 35,
                    "temp_outlet": 34,
                    "temp_source": 0,
                    "temp_test": -7
                },
                {"test_letter": "A",
                                                    "capacity": 8.8,
                                                    "cop": 3.2,
                                                    "design_flow_temp": 55,
                                                    "temp_outlet": 52,
                                                    "temp_source": 0,
                                                    "temp_test": -7
                                                   },
                                                   {
                                                    "test_letter": "B",
                                                    "capacity": 8.6,
                                                    "cop": 3.6,
                                                    "design_flow_temp": 55,
                                                    "temp_outlet": 42,
                                                    "temp_source": 0,
                                                    "temp_test": 2
                                                   },
                                                   {
                                                    "test_letter": "C",
                                                    "capacity": 8.5,
                                                    "cop": 3.9,
                                                    "design_flow_temp": 55,
                                                    "temp_outlet": 36,
                                                    "temp_source": 0,
                                                    "temp_test": 7
                                                   },
                                                   {
                                                    "test_letter": "D",
                                                    "capacity": 8.5,
                                                    "cop": 4.3,
                                                    "design_flow_temp": 55,
                                                    "temp_outlet": 30,
                                                    "temp_source": 0,
                                                    "temp_test": 12
                                                   },
                                                   {
                                                    "test_letter": "F",
                                                    "capacity": 8.8,
                                                    "cop": 3.2,
                                                    "design_flow_temp": 55,
                                                    "temp_outlet": 52,
                                                    "temp_source": 0,
                                                    "temp_test": -7
                                                   }
            ]
        });

        serde_json::from_value(input).unwrap()
    }

    fn create_heat_pump_with_exhaust_input_from_json(
        time_delay_backup: f64,
        temp_return_feed_max: u32,
    ) -> HeatSourceWetDetails {
        let input = json!({
            "type": "HeatPump",
            "EnergySupply": "mains gas",
            "source_type": "ExhaustAirMixed",
            "sink_type": "Water",
            "backup_ctrl_type": "Substitute",
            "temp_distribution_heat_network": 1.,
            "time_delay_backup": time_delay_backup,
            "modulating_control": true,
            "min_modulation_rate_20": 20.,
            "min_modulation_rate_35": 0.35,
            "min_modulation_rate_55": 0.4,
            "time_constant_onoff_operation": 140,
            "temp_return_feed_max": temp_return_feed_max,
            "temp_lower_operating_limit": -5.0,
            "min_temp_diff_flow_return_for_hp_to_operate": 0.0,
            "var_flow_temp_ctrl_during_test": true,
            "power_heating_circ_pump": 0.015,
            "power_source_circ_pump": 0.01,
            "power_standby": 0.015,
            "power_crankcase_heater": 0.01,
            "power_off": 0.015,
            "power_max_backup": 3.0,
            "eahp_mixed_max_temp": 10,
            "eahp_mixed_min_temp": 0,
            "test_data_EN14825": [
                {
                    "air_flow_rate": 100.0,
                    "test_letter": "A",
                    "capacity": 8.4,
                    "cop": 4.6,
                    "design_flow_temp": 35,
                    "temp_outlet": 34,
                    "temp_source": 0,
                    "temp_test": -7,
                    "eahp_mixed_ext_air_ratio": 0.62
                },
                {
                    "air_flow_rate": 100.0,
                    "test_letter": "B",
                    "capacity": 8.3,
                    "cop": 4.9,
                    "design_flow_temp": 35,
                    "temp_outlet": 30,
                    "temp_source": 0,
                    "temp_test": 2,
                    "eahp_mixed_ext_air_ratio": 0.62
                },
                {
                    "air_flow_rate": 100.0,
                    "test_letter": "C",
                    "capacity": 8.3,
                    "cop": 5.1,
                    "design_flow_temp": 35,
                    "temp_outlet": 27,
                    "temp_source": 0,
                    "temp_test": 7,
                    "eahp_mixed_ext_air_ratio": 0.62
                },
                {
                    "air_flow_rate": 100.0,
                    "test_letter": "D",
                    "capacity": 8.2,
                    "cop": 5.4,
                    "design_flow_temp": 35,
                    "temp_outlet": 24,
                    "temp_source": 0,
                    "temp_test": 12,
                    "eahp_mixed_ext_air_ratio": 0.62
                },
                {
                    "air_flow_rate": 100.0,
                    "test_letter": "F",
                    "capacity": 8.4,
                    "cop": 4.6,
                    "design_flow_temp": 35,
                    "temp_outlet": 34,
                    "temp_source": 0,
                    "temp_test": -7,
                    "eahp_mixed_ext_air_ratio": 0.62
                }
            ]
        });

        serde_json::from_value(input).unwrap()
    }

    fn create_heat_pump_sink_air_input_from_json() -> HeatSourceWetDetails {
        let input = json!({
            "type": "HeatPump",
            "EnergySupply": "mains gas",
            "source_type": "OutsideAir",
            "sink_type": "Air",
            "backup_ctrl_type": "TopUp",
            "time_delay_backup": 2.0,
            "modulating_control": true,
            "min_modulation_rate_35": 0.35,
            "min_modulation_rate_55": 0.4,
            "time_constant_onoff_operation": 140,
            "temp_return_feed_max": 70.0,
            "temp_lower_operating_limit": -5.0,
            "min_temp_diff_flow_return_for_hp_to_operate": 0.0,
            "var_flow_temp_ctrl_during_test": true,
            "power_heating_circ_pump": 0.015,
            "power_source_circ_pump": 0.01,
            "power_standby": 0.015,
            "power_crankcase_heater": 0.01,
            "power_off": 0.015,
            "power_max_backup": 3.0,
            "min_modulation_rate_20": 20.0,
            "test_data_EN14825": [
                {
                    "test_letter": "A",
                    "capacity": 8.4,
                    "cop": 4.6,
                    "design_flow_temp": 35,
                    "temp_outlet": 34,
                    "temp_source": 0,
                    "temp_test": -7
                },
                {
                    "test_letter": "B",
                    "capacity": 8.3,
                    "cop": 4.9,
                    "design_flow_temp": 35,
                    "temp_outlet": 30,
                    "temp_source": 0,
                    "temp_test": 2
                },
                {
                    "test_letter": "C",
                    "capacity": 8.3,
                    "cop": 5.1,
                    "design_flow_temp": 35,
                    "temp_outlet": 27,
                    "temp_source": 0,
                    "temp_test": 7
                },
                {
                    "test_letter": "D",
                    "capacity": 8.2,
                    "cop": 5.4,
                    "design_flow_temp": 35,
                    "temp_outlet": 24,
                    "temp_source": 0,
                    "temp_test": 12
                },
                {
                    "test_letter": "F",
                    "capacity": 8.4,
                    "cop": 4.6,
                    "design_flow_temp": 35,
                    "temp_outlet": 34,
                    "temp_source": 0,
                    "temp_test": -7
                }
            ]
        });
        serde_json::from_value(input).unwrap()
    }

    fn create_heat_pump_nw_input_from_json() -> HeatSourceWetDetails {
        let input = json!({
            "type": "HeatPump",
            "EnergySupply": "mains gas",
            "source_type": "HeatNetwork",
            "sink_type": "Water",
            "backup_ctrl_type": "TopUp",
            "temp_distribution_heat_network": 20.0,
            "time_delay_backup": 2.0,
            "modulating_control": true,
            "min_modulation_rate_35": 0.35,
            "min_modulation_rate_55": 0.4,
            "time_constant_onoff_operation": 140,
            "temp_return_feed_max": 70.0,
            "temp_lower_operating_limit": -5.0,
            "min_temp_diff_flow_return_for_hp_to_operate": 0.0,
            "var_flow_temp_ctrl_during_test": false,
            "power_heating_circ_pump": 0.015,
            "power_source_circ_pump": 0.01,
            "power_standby": 0.015,
            "power_crankcase_heater": 0.01,
            "power_off": 0.015,
            "power_max_backup": 3.0,
            "test_data_EN14825": [
                {
                    "test_letter": "A",
                    "capacity": 8.4,
                    "cop": 4.6,
                    "design_flow_temp": 35,
                    "temp_outlet": 34,
                    "temp_source": 0,
                    "temp_test": -7
                },
                {
                    "test_letter": "B",
                    "capacity": 8.3,
                    "cop": 4.9,
                    "design_flow_temp": 35,
                    "temp_outlet": 30,
                    "temp_source": 0,
                    "temp_test": 2
                },
                {
                    "test_letter": "C",
                    "capacity": 8.3,
                    "cop": 5.1,
                    "design_flow_temp": 35,
                    "temp_outlet": 27,
                    "temp_source": 0,
                    "temp_test": 7
                },
                {
                    "test_letter": "D",
                    "capacity": 8.2,
                    "cop": 5.4,
                    "design_flow_temp": 35,
                    "temp_outlet": 24,
                    "temp_source": 0,
                    "temp_test": 12
                },
                {
                    "test_letter": "F",
                    "capacity": 8.4,
                    "cop": 4.6,
                    "design_flow_temp": 35,
                    "temp_outlet": 34,
                    "temp_source": 0,
                    "temp_test": -7
                }
            ]
        });
        serde_json::from_value(input).unwrap()
    }

    fn create_heat_pump(
        heat_pump_input: HeatSourceWetDetails,
        energy_supply_conn_name_auxiliary: &str,
        energy_supply_heat_source: Option<Arc<RwLock<EnergySupply>>>,
        boiler: Option<Arc<RwLock<Boiler>>>,
        throughput_exhaust_air: Option<f64>,
        temp_internal_air: Option<f64>,
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
        cost_schedule_hybrid_hp: Option<Value>,
    ) -> HeatPump {
        let energy_supply = energy_supply(simulation_time_for_heat_pump);
        let number_of_zones = 2;
        let cost_schedule_hybrid_hp =
            cost_schedule_hybrid_hp.map(|json| serde_json::from_value(json).unwrap());
        let temp_internal_air_fn = create_temp_internal_air_fn(temp_internal_air.unwrap_or(20.));
        HeatPump::new(
            &heat_pump_input,
            Arc::from(RwLock::from(energy_supply)),
            energy_supply_conn_name_auxiliary,
            simulation_time_for_heat_pump.step,
            Arc::new(external_conditions),
            number_of_zones,
            throughput_exhaust_air,
            energy_supply_heat_source,
            false,
            boiler,
            cost_schedule_hybrid_hp,
            temp_internal_air_fn,
        )
        .unwrap()
    }

    fn create_default_heat_pump(
        energy_supply_conn_name_auxiliary: Option<&str>,
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
        backup_ctrl_type: Option<&str>,
    ) -> HeatPump {
        let energy_supply_conn_name_auxiliary =
            energy_supply_conn_name_auxiliary.unwrap_or("HeatPump_auxiliary: hp");

        let heat_pump_input = create_heat_pump_input_from_json(backup_ctrl_type, None);

        create_heat_pump(
            heat_pump_input,
            energy_supply_conn_name_auxiliary,
            None, // energy_supply_heat_source: Option<Arc<RwLock<EnergySupply>>>,
            None, // boiler: Option<Arc<Mutex<Boiler>>>,
            None, // throughput_exhaust_air: Option<f64>,
            None, // temp_internal_air: Option<f64>,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        )
    }

    fn create_heat_pump_with_boiler(
        energy_supply_conn_name_auxiliary: &str,
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) -> HeatPump {
        let boiler = Arc::from(RwLock::from(create_boiler(
            external_conditions.clone(),
            energy_supply(simulation_time_for_heat_pump),
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));

        let input = create_heat_pump_input_from_json(None, None);

        create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            Some(boiler),
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        )
    }

    fn create_heat_pump_with_exhaust(
        energy_supply_conn_name_auxiliary: &str,
        time_delay_backup: Option<f64>,
        temp_return_feed_max: Option<u32>,
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) -> HeatPump {
        let time_delay_backup = time_delay_backup.unwrap_or(2.);
        let temp_return_feed_max = temp_return_feed_max.unwrap_or(70);
        let heat_pump_input =
            create_heat_pump_with_exhaust_input_from_json(time_delay_backup, temp_return_feed_max);

        create_heat_pump(
            heat_pump_input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            Some(101.), // throughput_exhaust_air
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        )
    }

    fn create_heat_pump_sink_air(
        energy_supply_conn_name_auxiliary: &str,
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) -> HeatPump {
        let input = create_heat_pump_sink_air_input_from_json();
        create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        )
    }

    fn create_heat_pump_with_energy_supply_heat_source(
        energy_supply_conn_name_auxiliary: &str,
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) -> HeatPump {
        let energy_supply_heat_source = RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Custom,
                simulation_time_for_heat_pump.iter().total_steps(),
            )
            .build(),
        );

        let input = create_heat_pump_nw_input_from_json();
        create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            Some(energy_supply_heat_source.into()),
            None,
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        )
    }

    fn create_setpoint_time_control(schedule: Vec<Option<f64>>) -> Control {
        Control::SetpointTime(SetpointTimeControl::new(
            schedule,
            0,
            1.,
            Default::default(),
            Default::default(),
            1.,
        ))
    }

    #[rstest]
    fn test_init_no_backup_ctrl_type(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let input = create_heat_pump_input_from_json(Some("None"), None);
        let heat_pump = create_heat_pump(
            input,
            "HeatPump_auxiliary: hp",
            None,
            None,
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        assert_eq!(heat_pump.backup_heater_max_power(), 0.);
    }

    // skipping Python's test_init_missing_heat_network as the Rust would error earlier and not reach the constructor/new function

    #[rstest]
    #[should_panic]
    fn test_init_overvent_ratio_too_high(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        // Test that the constructor throws if the overvent ratio is above 1
        let mut input = create_heat_pump_input_from_json(None, Some("ExhaustAirMVHR"));
        if let HeatSourceWetDetails::HeatPump {
            ref mut test_data_en14825,
            ..
        } = input
        {
            *test_data_en14825 = vec![
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(1050.),
                    test_letter: TestLetter::A,
                    capacity: 10.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 40.,
                    temp_test: 40.,
                    eahp_mixed_ext_air_ratio: None,
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(1050.),
                    test_letter: TestLetter::B,
                    capacity: 11.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 41.,
                    temp_test: 41.,
                    eahp_mixed_ext_air_ratio: None,
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(1050.),
                    test_letter: TestLetter::C,
                    capacity: 12.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 42.,
                    temp_test: 42.,
                    eahp_mixed_ext_air_ratio: None,
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(1050.),
                    test_letter: TestLetter::D,
                    capacity: 13.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 43.,
                    temp_test: 43.,
                    eahp_mixed_ext_air_ratio: None,
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(1050.),
                    test_letter: TestLetter::F,
                    capacity: 14.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 44.,
                    temp_test: 44.,
                    eahp_mixed_ext_air_ratio: None,
                },
            ];
        }
        create_heat_pump(
            input,
            "HeatPump_auxiliary: hp",
            None,
            None,
            Some(1000.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
    }

    #[rstest]
    #[should_panic] // TODO this errors for a different reason to the Python and seems to relate to this comment above:
                    // in the Python there is a check here about an air flow rate being present - decided
                    // any validation should really happen upstream of here
    fn test_init_air_flow_rate(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut input = create_heat_pump_input_from_json(None, Some("Ground"));
        if let HeatSourceWetDetails::HeatPump {
            ref mut test_data_en14825,
            ..
        } = input
        {
            *test_data_en14825 = vec![
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::A,
                    capacity: 14.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 44.,
                    temp_test: 44.,
                    eahp_mixed_ext_air_ratio: None,
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::B,
                    capacity: 14.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 44.,
                    temp_test: 44.,
                    eahp_mixed_ext_air_ratio: None,
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::C,
                    capacity: 14.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 44.,
                    temp_test: 44.,
                    eahp_mixed_ext_air_ratio: None,
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::D,
                    capacity: 14.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 44.,
                    temp_test: 44.,
                    eahp_mixed_ext_air_ratio: None,
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::F,
                    capacity: 14.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 44.,
                    temp_test: 44.,
                    eahp_mixed_ext_air_ratio: None,
                },
            ];
        }
        create_heat_pump(
            input,
            "HeatPump_auxiliary: hp",
            None,
            None,
            Some(1000.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
    }

    #[rstest]
    #[should_panic]
    fn test_init_unique_air_ratio(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut input = create_heat_pump_input_from_json(None, Some("ExhaustAirMixed"));
        if let HeatSourceWetDetails::HeatPump {
            ref mut eahp_mixed_min_temp,
            ref mut eahp_mixed_max_temp,
            ref mut test_data_en14825,
            ..
        } = input
        {
            *eahp_mixed_min_temp = Some(20.);
            *eahp_mixed_max_temp = Some(30.);
            *test_data_en14825 = vec![
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::A,
                    capacity: 10.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 40.,
                    temp_test: 40.,
                    eahp_mixed_ext_air_ratio: Some(0.5),
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::B,
                    capacity: 11.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 41.,
                    temp_test: 41.,
                    eahp_mixed_ext_air_ratio: Some(0.5),
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::C,
                    capacity: 12.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 42.,
                    temp_test: 42.,
                    eahp_mixed_ext_air_ratio: Some(0.5),
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::D,
                    capacity: 13.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 43.,
                    temp_test: 43.,
                    eahp_mixed_ext_air_ratio: Some(0.5),
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::F,
                    capacity: 14.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 44.,
                    temp_test: 44.,
                    eahp_mixed_ext_air_ratio: Some(0.6),
                },
            ];
        }
        create_heat_pump(
            input,
            "HeatPump_auxiliary: hp",
            None,
            None,
            Some(1000.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
    }

    #[rstest]
    fn test_source_is_exhaust_air(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
        heat_pump.source_type = HeatPumpSourceType::OutsideAir;

        assert!(!heat_pump.source_is_exhaust_air());
    }
    #[rstest]
    fn test_init_buffer_tank(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hp";
        let mut input = create_heat_pump_input_from_json(None, None);
        if let HeatSourceWetDetails::HeatPump {
            ref mut buffer_tank,
            ..
        } = input
        {
            *buffer_tank = Some(Box::new(HeatPumpBufferTank {
                daily_losses: 10.,
                volume: 100.,
                pump_fixed_flow_rate: 5.,
                pump_power_at_flow_rate: 10.,
            }));
        }

        let heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            Some(1000.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        assert!(heat_pump.buffer_tank.is_some());
        assert_eq!(heat_pump.buffer_tank.unwrap().volume, 100.);
    }

    #[rstest]
    fn test_buffer_int_gains(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hp";
        let mut input = create_heat_pump_input_from_json(None, None);
        if let HeatSourceWetDetails::HeatPump {
            ref mut buffer_tank,
            ..
        } = input
        {
            *buffer_tank = Some(Box::new(HeatPumpBufferTank {
                daily_losses: 10.,
                volume: 100.,
                pump_fixed_flow_rate: 5.,
                pump_power_at_flow_rate: 10.,
            }));
        }

        let mut heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            Some(1000.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
        if let Some(buffer_tank) = heat_pump.buffer_tank.as_mut() {
            buffer_tank.q_heat_loss_buffer_rbl = 0.0015; // This makes buffer_tank.internal_gains() return 3.
                                                         // Python does this: self.heat_pump._HeatPump__buffer_tank.internal_gains.return_value = 3
            assert_eq!(buffer_tank.internal_gains(), 3.);
        }

        assert_eq!(heat_pump.buffer_int_gains(), 3.);
    }

    #[rstest]
    fn test_buffer_int_gains_no_buffer_tank(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hp";
        let input = create_heat_pump_input_from_json(None, None);
        let heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            Some(1000.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        assert_eq!(heat_pump.buffer_int_gains(), 0.);
    }

    #[rstest]
    fn test_create_service_connection(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let service_name = "new_service";
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );
        let heat_pump = Arc::from(Mutex::from(heat_pump));

        // Ensure the service name does not exist in energy_supply_connections
        assert!(!heat_pump
            .lock()
            .energy_supply_connections
            .contains_key(service_name));

        // Call the method under test
        HeatPump::create_service_connection(heat_pump.clone(), service_name).unwrap();

        // Check that the service name was added to energy_supply_connections
        assert!(heat_pump
            .lock()
            .energy_supply_connections
            .contains_key(service_name));

        let create_connection_result = HeatPump::create_service_connection(heat_pump, service_name);

        // second attempt to create a service connection with same name should error
        assert!(create_connection_result.is_err());

        // Check with heat_network
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hp1";

        let heat_pump_nw = create_heat_pump_with_energy_supply_heat_source(
            energy_supply_conn_name_auxiliary,
            external_conditions,
            simulation_time_for_heat_pump,
        );
        let heat_pump_nw = Arc::from(Mutex::from(heat_pump_nw));

        HeatPump::create_service_connection(heat_pump_nw.clone(), "new_service_nw").unwrap();

        assert!(heat_pump_nw
            .lock()
            .energy_supply_heat_source_connections
            .contains_key("new_service_nw"));
    }

    #[rstest]
    fn test_create_service_hot_water_combi(
        simulation_time_for_heat_pump: SimulationTime,
        external_conditions: ExternalConditions,
    ) {
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: boiler";
        let heat_pump_with_boiler = create_heat_pump_with_boiler(
            energy_supply_conn_name_auxiliary,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
        );
        let boiler_data: HotWaterSourceDetails = HotWaterSourceDetails::CombiBoiler {
            cold_water_source: "mains water".into(),
            heat_source_wet: "hp".into(),
            separate_dhw_tests: BoilerHotWaterTest::ML,
            rejected_energy_1: Some(0.0004),
            storage_loss_factor_2: Some(0.91574),
            rejected_factor_3: Some(0.),
            setpoint_temp: Some(60.),
            daily_hw_usage: 120.,
        };
        let service_name = "service_hot_water_combi";
        let temp_hot_water = 50.;
        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(
            ColdWaterSource::new(vec![1.0, 1.2], 0, simulation_time_for_heat_pump.step).into(),
        );
        let boiler_service_water_combi: Result<BoilerServiceWaterCombi, anyhow::Error> =
            heat_pump_with_boiler.create_service_hot_water_combi(
                boiler_data.clone(),
                service_name,
                temp_hot_water,
                cold_feed.clone(),
            );

        assert!(matches!(
            boiler_service_water_combi,
            Ok(BoilerServiceWaterCombi { .. })
        ));

        let heat_pump = create_default_heat_pump(
            Some("HeatPump_auxiliary: boiler1"),
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        let boiler_service_water_combi: Result<BoilerServiceWaterCombi, anyhow::Error> = heat_pump
            .create_service_hot_water_combi(boiler_data, service_name, temp_hot_water, cold_feed);

        // creating a BoilerServiceWaterCombi should error on heat pump without boiler
        assert!(boiler_service_water_combi.is_err());
    }

    #[rstest]
    fn test_create_service_hot_water(
        simulation_time_for_heat_pump: SimulationTime,
        external_conditions: ExternalConditions,
    ) {
        let heat_pump = create_default_heat_pump(
            Some("HeatPump_auxiliary: HotWater"),
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let heat_pump = Arc::from(Mutex::from(heat_pump));

        let service_name = "service_hot_water";

        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(
            ColdWaterSource::new(vec![1.0, 1.2], 0, simulation_time_for_heat_pump.step).into(),
        );

        let control_min = Arc::new(create_setpoint_time_control(vec![
            Some(52.),
            Some(52.),
            None,
            Some(52.),
        ]));
        let control_max = Arc::new(create_setpoint_time_control(vec![
            Some(60.),
            Some(60.),
            Some(60.),
            Some(60.),
        ]));

        let hot_water_service = HeatPump::create_service_hot_water(
            heat_pump.clone(),
            service_name,
            60.,
            cold_feed.clone().into(),
            control_min.clone(),
            control_max.clone(),
        )
        .unwrap();

        assert!(matches!(hot_water_service, HeatPumpServiceWater { .. }));

        assert!(heat_pump
            .lock()
            .energy_supply_connections
            .contains_key(service_name));

        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: HotWater_boiler";

        let heat_pump_with_boiler = Mutex::from(create_heat_pump_with_boiler(
            energy_supply_conn_name_auxiliary,
            external_conditions,
            simulation_time_for_heat_pump,
        ));

        let hot_water_service = HeatPump::create_service_hot_water(
            heat_pump_with_boiler.into(),
            service_name,
            60.,
            cold_feed.into(),
            control_min,
            control_max,
        )
        .unwrap();

        assert!(hot_water_service.hybrid_boiler_service.is_some());
    }

    #[rstest]
    fn test_create_service_space_heating(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let service_name = "service_space";
        let temp_limit_upper = 50.0;
        let temp_diff_emit_design = 50.0;
        let control = Arc::from(Control::OnOffTime(OnOffTimeControl::new(vec![], 0, 0.)));
        let volume_heated = 250.0;

        let heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let heat_pump = Arc::from(Mutex::from(heat_pump));

        let service_space_heating = HeatPump::create_service_space_heating(
            heat_pump.clone(),
            service_name,
            HeatPumpEmitterType::RadiatorsUfh,
            temp_limit_upper,
            temp_diff_emit_design,
            55.,
            control.clone(),
            volume_heated,
        );

        assert!(matches!(service_space_heating, HeatPumpServiceSpace { .. }));

        assert!(heat_pump
            .lock()
            .energy_supply_connections
            .contains_key(service_name));

        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: SpaceHeating_boiler";
        let heat_pump_with_boiler = create_heat_pump_with_boiler(
            energy_supply_conn_name_auxiliary,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
        );

        let heat_pump_with_boiler = Arc::from(Mutex::from(heat_pump_with_boiler));

        let service_name = "service_space_boiler";

        HeatPump::create_service_space_heating(
            heat_pump_with_boiler.clone(),
            service_name,
            HeatPumpEmitterType::RadiatorsUfh,
            temp_limit_upper,
            temp_diff_emit_design,
            55.,
            control.clone(),
            volume_heated,
        );

        assert!(heat_pump_with_boiler
            .lock()
            .energy_supply_connections
            .contains_key(service_name));

        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: exhaust_service_space";

        let heat_pump_exhaust = create_heat_pump_with_exhaust(
            energy_supply_conn_name_auxiliary,
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
        );

        let heat_pump_exhaust = Arc::from(Mutex::from(heat_pump_exhaust));

        let service_name = "service_space_exhaust";

        HeatPump::create_service_space_heating(
            heat_pump_exhaust.clone(),
            service_name,
            HeatPumpEmitterType::RadiatorsUfh,
            temp_limit_upper,
            temp_diff_emit_design,
            55.,
            control,
            volume_heated,
        );

        assert_relative_eq!(
            heat_pump_exhaust.lock().volume_heated_all_services.unwrap(),
            250.
        );
    }

    #[rstest]
    fn test_create_service_space_heating_warm_air(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );
        let heat_pump = Arc::from(Mutex::from(heat_pump));

        let service_name = "service_space_warmair";
        let control = Arc::from(Control::OnOffTime(OnOffTimeControl::new(vec![], 0, 0.)));
        let volume_heated = 250.;
        let frac_convective = 0.9;

        let result = HeatPump::create_service_space_heating_warm_air(
            heat_pump,
            service_name,
            control.clone(),
            frac_convective,
            volume_heated,
        );

        // Check we get an error when Sink type is not air
        assert!(result.is_err());

        // Check without boiler object and sink type 'AIR'
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: sink_air";
        let heat_pump_sink_air = create_heat_pump_sink_air(
            energy_supply_conn_name_auxiliary,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
        );
        let heat_pump_sink_air = Arc::from(Mutex::from(heat_pump_sink_air));

        let result = HeatPump::create_service_space_heating_warm_air(
            heat_pump_sink_air.clone(),
            service_name,
            control.clone(),
            frac_convective,
            volume_heated,
        );

        assert!(result.is_ok());
        assert!(heat_pump_sink_air
            .lock()
            .energy_supply_connections
            .contains_key(service_name));

        // Check that warm air hybrid heat pump service cannot be created
        let boiler = Arc::from(RwLock::from(create_boiler(
            external_conditions.clone(),
            energy_supply(simulation_time_for_heat_pump),
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));

        let input = create_heat_pump_sink_air_input_from_json();
        let heat_pump_sink_air_with_boiler = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            Some(boiler),
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
        let result = HeatPump::create_service_space_heating_warm_air(
            Arc::new(Mutex::new(heat_pump_sink_air_with_boiler)),
            service_name,
            control.clone(),
            frac_convective,
            volume_heated,
        );

        assert!(result.is_err());
    }

    #[rstest]
    fn test_create_service_space_heating_warm_air_exhaust_air(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut input = create_heat_pump_input_from_json(None, Some("ExhaustAirMixed"));
        if let HeatSourceWetDetails::HeatPump {
            ref mut eahp_mixed_min_temp,
            ref mut eahp_mixed_max_temp,
            ref mut test_data_en14825,
            ..
        } = input
        {
            *eahp_mixed_min_temp = Some(20.);
            *eahp_mixed_max_temp = Some(30.);
            *test_data_en14825 = vec![
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::A,
                    capacity: 10.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 40.,
                    temp_test: 40.,
                    eahp_mixed_ext_air_ratio: Some(0.5),
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::B,
                    capacity: 11.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 41.,
                    temp_test: 41.,
                    eahp_mixed_ext_air_ratio: Some(0.5),
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::C,
                    capacity: 12.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 42.,
                    temp_test: 42.,
                    eahp_mixed_ext_air_ratio: Some(0.5),
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::D,
                    capacity: 13.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 43.,
                    temp_test: 43.,
                    eahp_mixed_ext_air_ratio: Some(0.5),
                },
                input::HeatPumpTestDatum {
                    air_flow_rate: Some(50.),
                    test_letter: TestLetter::F,
                    capacity: 14.,
                    cop: 1.,
                    design_flow_temp: 45.,
                    temp_outlet: 35.,
                    temp_source: 44.,
                    temp_test: 44.,
                    eahp_mixed_ext_air_ratio: Some(0.5),
                },
            ];
        }
        let mut heat_pump = create_heat_pump(
            input,
            "HeatPump_auxiliary: hp",
            None,
            None,
            Some(1000.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        heat_pump.sink_type = HeatPumpSinkType::Air;
        let heat_pump = Arc::from(Mutex::from(heat_pump));

        let control = Arc::from(Control::OnOffTime(OnOffTimeControl::new(vec![], 0, 0.)));
        // in Python None is passed in for the control but this is not optional in Rust
        HeatPump::create_service_space_heating_warm_air(
            heat_pump.clone(),
            "mains_gas",
            control,
            0.5,
            100.,
        )
        .unwrap();

        assert_eq!(heat_pump.lock().volume_heated_all_services, Some(100.));
    }

    #[rstest]
    fn test_get_temp_source(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let result =
            heat_pump.get_temp_source(simulation_time_for_heat_pump.iter().current_iteration());

        assert_relative_eq!(result, 273.15);

        // Check with ExhaustAirMixed
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: exhaust_source";
        let heat_pump_with_exhaust = create_heat_pump_with_exhaust(
            energy_supply_conn_name_auxiliary,
            None,
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
        );

        let result = heat_pump_with_exhaust
            .get_temp_source(simulation_time_for_heat_pump.iter().current_iteration());

        assert_relative_eq!(result, 280.75);

        // Check with heat_network
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: heat_nw";

        let heat_pump_nw = create_heat_pump_with_energy_supply_heat_source(
            energy_supply_conn_name_auxiliary,
            external_conditions,
            simulation_time_for_heat_pump,
        );

        let result =
            heat_pump_nw.get_temp_source(simulation_time_for_heat_pump.iter().current_iteration());

        assert_relative_eq!(result, 293.15);
    }

    #[rstest]
    fn test_get_temp_source_waterground(
        mut external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        external_conditions.air_temps = vec![25.; 8760]; // Equivalent of Python: self.extcond.air_temp_annual.return_value = 25
        let input = create_heat_pump_input_from_json(None, Some("WaterGround"));
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: heat_nw";

        let heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            None,
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        assert_relative_eq!(
            heat_pump.get_temp_source(simulation_time_for_heat_pump.iter().current_iteration()),
            298.15
        );
    }

    #[rstest]
    fn test_get_temp_source_watersurface(
        mut external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        external_conditions.air_temps = vec![26.; 8760]; // Equivalent of Python: self.extcond.air_temp_monthly.return_value = 26
        let input = create_heat_pump_input_from_json(None, Some("WaterSurface"));
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: heat_nw";

        let heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        assert_relative_eq!(
            heat_pump.get_temp_source(simulation_time_for_heat_pump.iter().current_iteration()),
            299.15
        );
    }

    #[rstest]
    fn test_get_temp_source_airmvhr(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let input = create_heat_pump_input_from_json(None, Some("ExhaustAirMVHR"));
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: heat_nw";

        let heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            None,
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        assert_relative_eq!(
            heat_pump.get_temp_source(simulation_time_for_heat_pump.iter().current_iteration()),
            293.15
        );
    }

    #[rstest]
    fn test_get_temp_source_ground(
        mut external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        external_conditions.air_temps = vec![27.; 8760]; // Equivalent of Python: self.extcond.air_temp.return_value = 27
        let input = create_heat_pump_input_from_json(None, Some("Ground"));
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: heat_nw";

        let heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        assert_relative_eq!(
            heat_pump.get_temp_source(simulation_time_for_heat_pump.iter().current_iteration()),
            281.15
        );
    }

    #[rstest]
    fn test_get_temp_source_airmev(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let input = create_heat_pump_input_from_json(None, Some("ExhaustAirMEV"));
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: heat_nw";

        let heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            None,
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        assert_relative_eq!(
            heat_pump.get_temp_source(simulation_time_for_heat_pump.iter().current_iteration()),
            293.15
        );
    }

    #[rstest]
    fn test_thermal_capacity_op_cond(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        // Check with source_type OutsideAir
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );
        let result = heat_pump
            .thermal_capacity_op_cond(290., 260., design_flow_temp_op_cond_k(55.))
            .unwrap();
        assert_relative_eq!(result, 10.69686174950112);

        // Check with ExhaustAirMixed
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: exhaust_source_capacity";
        let heat_pump_with_exhaust = create_heat_pump_with_exhaust(
            energy_supply_conn_name_auxiliary,
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
        );
        let result = heat_pump_with_exhaust
            .thermal_capacity_op_cond(300., 270., design_flow_temp_op_cond_k(55.))
            .unwrap();
        assert_relative_eq!(result, 8.70672362099314);
    }

    #[rstest]
    fn test_backup_energy_output_max(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        // With TopUp backup control
        let temp_output = 300.;
        let temp_return_feed = 290.;
        let time_available = 1.;
        let time_start = 0.;

        let heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let result = heat_pump
            .backup_energy_output_max(
                temp_output,
                temp_return_feed,
                time_available,
                time_start,
                None,
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        assert_relative_eq!(result, 3.);

        // With substitute backup control
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: exhaust_backup_energy";

        let heat_pump_with_exhaust = create_heat_pump_with_exhaust(
            energy_supply_conn_name_auxiliary,
            None,
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
        );

        let result = heat_pump_with_exhaust
            .backup_energy_output_max(
                temp_output,
                temp_return_feed,
                time_available,
                time_start,
                None,
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        assert_relative_eq!(result, 3.);

        // With a hybrid boiler
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hybrid_boiler";

        let boiler = Arc::from(RwLock::from(create_boiler(
            external_conditions.clone(),
            energy_supply(simulation_time_for_heat_pump),
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));

        let control = Arc::from(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(21.), Some(22.)],
            0,
            1.,
            Default::default(),
            Default::default(),
            simulation_time_for_heat_pump.step,
        )));

        let boiler_service_space = Boiler::create_service_space_heating(
            boiler.clone(),
            "service_boilerspace",
            control.clone(),
        );
        let hybrid_boiler_service =
            HybridBoilerService::Space(Arc::from(Mutex::from(boiler_service_space)));
        let input = create_heat_pump_input_from_json(None, None);

        let heat_pump_with_boiler = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            Some(boiler.clone()),
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        let result = heat_pump_with_boiler
            .backup_energy_output_max(
                temp_output,
                temp_return_feed,
                time_available,
                time_start,
                Some(hybrid_boiler_service),
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        assert_relative_eq!(result, 24.);

        // Test with boiler service water regular
        let boiler_service_water_regular = Boiler::create_service_hot_water_regular(
            boiler,
            "service_boilerwater",
            control.clone(),
            control,
        )
        .unwrap();
        let hybrid_boiler_service =
            HybridBoilerService::Regular(Arc::from(Mutex::from(boiler_service_water_regular)));

        assert_relative_eq!(
            heat_pump
                .backup_energy_output_max(
                    temp_output,
                    temp_return_feed,
                    time_available,
                    time_start,
                    Some(hybrid_boiler_service),
                    simulation_time_for_heat_pump.iter().current_iteration()
                )
                .unwrap(),
            24.
        );
    }

    #[rstest]
    fn test_energy_output_max(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let temp_output = 300.;
        let temp_return_feed = 295.;
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hp";
        let mut heat_pump = create_default_heat_pump(
            Some(energy_supply_conn_name_auxiliary),
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );
        let boiler = Arc::from(RwLock::from(create_boiler(
            external_conditions.clone(),
            energy_supply(simulation_time_for_heat_pump),
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));

        let control = Arc::from(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(0.), Some(0.)],
            0,
            1.,
            Default::default(),
            Default::default(),
            simulation_time_for_heat_pump.step,
        )));

        let hybrid_boiler_service = HybridBoilerService::Space(Arc::new(Mutex::new(
            BoilerServiceSpace::new(boiler.clone(), "boiler_service_space".into(), control),
        )));

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            assert_relative_eq!(
                heat_pump
                    .energy_output_max(
                        temp_output,
                        temp_return_feed,
                        design_flow_temp_op_cond_k(55.),
                        Some(hybrid_boiler_service.clone()),
                        None,
                        None,
                        None,
                        None,
                        Some("hp_service"),
                        t_it
                    )
                    .unwrap()
                    .0,
                [11.203924743692589, 11.51438002341549][t_idx]
            );
        }

        // Check with backup SUBSTITUTE
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: backup_substitute";
        let input = create_heat_pump_with_exhaust_input_from_json(2., 70);
        let mut heat_pump_exhaust = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            Some(101.),
            Some(30.),
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            assert_relative_eq!(
                heat_pump_exhaust
                    .energy_output_max(
                        temp_output,
                        temp_return_feed,
                        design_flow_temp_op_cond_k(55.),
                        Some(hybrid_boiler_service.clone()),
                        None,
                        None,
                        None,
                        None,
                        Some("heat_pump_exhaust_service"),
                        t_it
                    )
                    .unwrap()
                    .0,
                [10.191526455038767, 10.358981076999134][t_idx]
            );
        }
    }

    #[rstest]
    fn test_energy_output_max_buffer_tank(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let temp_output = 300.;
        let temp_return_feed = 295.;
        let emitters_data_for_buffer_tank = EMITTERS_DATA_FOR_BUFFER_TANK[0];
        let mut input = create_heat_pump_input_from_json(None, None);
        if let HeatSourceWetDetails::HeatPump {
            ref mut buffer_tank,
            ..
        } = input
        {
            *buffer_tank = Some(Box::new(HeatPumpBufferTank {
                daily_losses: 10.,
                volume: 100.,
                pump_fixed_flow_rate: 50.,
                pump_power_at_flow_rate: 10.,
            }));
        }

        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hp";

        let mut heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            Some(1000.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            let (energy, emitters_data_for_buffer_tank) = heat_pump
                .energy_output_max(
                    temp_output,
                    temp_return_feed,
                    emitters_data_for_buffer_tank.design_flow_temp,
                    None,
                    None,
                    None,
                    None,
                    Some(emitters_data_for_buffer_tank),
                    Some("boiler_service_space"),
                    t_it,
                )
                .unwrap();

            assert_relative_eq!(energy, [8.2191174037396, 8.449538231375868][t_idx]);

            let result = emitters_data_for_buffer_tank.unwrap().result;

            assert_eq!(result.flow_temp_increase_due_to_buffer, 8.185825328614797);
            assert_eq!(result.pump_power_at_flow_rate, 10.);
            assert_eq!(result.heat_loss_buffer_kwh, 0.0964687074719922);
        }
    }

    #[rstest]
    fn test_energy_output_max_cost_schedule_hybrid_hp(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let temp_output = 300.;
        let temp_return_feed = 295.;
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hybrid_boiler";

        let boiler = Arc::from(RwLock::from(create_boiler(
            external_conditions.clone(),
            energy_supply(simulation_time_for_heat_pump),
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));

        let control = Arc::from(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(0.), Some(0.)],
            0,
            1.,
            Default::default(),
            Default::default(),
            simulation_time_for_heat_pump.step,
        )));

        let hybrid_boiler_service = HybridBoilerService::Space(Arc::new(Mutex::new(
            BoilerServiceSpace::new(boiler.clone(), "boiler_service_space".into(), control),
        )));

        let cost_schedule_hybrid_hp = json!({
           "cost_schedule_start_day": 0,
           "cost_schedule_time_series_step": 1,
           "cost_schedule_hp":
               {"main": [16, 18, 24, {"value": 16, "repeat": 5}]},
           "cost_schedule_boiler": {"main": [{"value": 4, "repeat": 8}]}
        });

        let input = create_heat_pump_input_from_json(None, None);

        let mut heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            Some(boiler),
            Some(1000.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            Some(cost_schedule_hybrid_hp),
        );

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            let energy = heat_pump
                .energy_output_max(
                    temp_output,
                    temp_return_feed,
                    design_flow_temp_op_cond_k(55.),
                    Some(hybrid_boiler_service.clone()),
                    None,
                    None,
                    None,
                    None,
                    Some("boiler_service_space"),
                    t_it,
                )
                .unwrap()
                .0;

            assert_relative_eq!(energy, [11.203924743692589, 11.51438002341549][t_idx]);
        }
    }

    #[rstest]
    fn test_energy_output_max_cost_schedule_hybrid_hp_not_cost_effective(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let temp_output = 300.;
        let temp_return_feed = 295.;
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hybrid_boiler";

        let boiler = Arc::from(RwLock::from(create_boiler(
            external_conditions.clone(),
            energy_supply(simulation_time_for_heat_pump),
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));

        let control = Arc::from(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(0.), Some(0.)],
            0,
            1.,
            Default::default(),
            Default::default(),
            simulation_time_for_heat_pump.step,
        )));

        let hybrid_boiler_service = HybridBoilerService::Space(Arc::new(Mutex::new(
            BoilerServiceSpace::new(boiler.clone(), "boiler_service_space".into(), control),
        )));

        let cost_schedule_hybrid_hp = json!({
           "cost_schedule_start_day": 0,
           "cost_schedule_time_series_step": 1,
           "cost_schedule_hp":
               {"main": [{"value": 40, "repeat": 8}]},
           "cost_schedule_boiler": {"main": [{"value": 5, "repeat": 8}]}
        });

        let input = create_heat_pump_input_from_json(None, None);

        let mut heat_pump = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            Some(boiler),
            Some(1000.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            Some(cost_schedule_hybrid_hp),
        );

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            let energy = heat_pump
                .energy_output_max(
                    temp_output,
                    temp_return_feed,
                    design_flow_temp_op_cond_k(55.),
                    Some(hybrid_boiler_service.clone()),
                    None,
                    None,
                    None,
                    None,
                    Some("boiler_service_space"),
                    t_it,
                )
                .unwrap()
                .0;

            assert_relative_eq!(energy, [24., 24.][t_idx]);
        }
    }

    // skipping Python's test_energy_output_max_topup due to mocking (self.extcond.air_temp.return_value = -5)

    #[rstest]
    fn test_cop_op_cond(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        // Check with type service type SPACE
        let temp_spread_correction = TempSpreadCorrectionArg::Float(1.);
        let service_type = HeatingServiceType::Space;
        let temp_output = 320.;
        let temp_source = 275.;

        let heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let cop_op_cond = heat_pump
            .cop_op_cond(
                Some(service_type),
                temp_output, // Kelvin
                temp_source, // Kelvin
                temp_spread_correction,
                design_flow_temp_op_cond_k(55.),
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        assert_relative_eq!(cop_op_cond, 3.4835975137798934);

        // Check with sink type 'AIR'
        let temp_spread_correction = TempSpreadCorrectionArg::Float(1.);
        let energy_supply_conn_name_auxiliary = "auxillary_cop_eff";

        let heat_pump_sink_air = create_heat_pump_sink_air(
            energy_supply_conn_name_auxiliary,
            external_conditions,
            simulation_time_for_heat_pump,
        );

        let cop_op_cond = heat_pump_sink_air
            .cop_op_cond(
                Some(service_type),
                temp_output,
                temp_source,
                temp_spread_correction,
                design_flow_temp_op_cond_k(55.),
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        assert_relative_eq!(cop_op_cond, 4.384741326956177);
    }

    #[rstest]
    fn test_cop_op_cond_lr_equal(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
        for value in heat_pump.test_data.test_data.values_mut().flatten() {
            value.theoretical_load_ratio = 1.1;
        }

        let temp_spread_correction = TempSpreadCorrectionArg::Float(1.);
        let service_type = HeatingServiceType::Space;
        let temp_output = 320.;
        let temp_source = 275.;

        let cop_op_cond = heat_pump
            .cop_op_cond(
                Some(service_type),
                temp_output,
                temp_source,
                temp_spread_correction,
                design_flow_temp_op_cond_k(55.),
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        assert_relative_eq!(cop_op_cond, 3.2295002426006794);
    }

    #[rstest]
    /// Check energy output limited by upper temperature
    fn test_energy_output_limited(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let energy_output_required = 1.5;
        let temp_output = 320.;
        let temp_used_for_scaling = 310.;
        let temp_limit_upper = 340.;

        let heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        let result = heat_pump.energy_output_limited(
            energy_output_required,
            Some(temp_output),
            temp_used_for_scaling,
            temp_limit_upper,
        );

        assert_relative_eq!(result, 1.5);

        let energy_output_required = 1.5;
        let temp_output = 320.;
        let temp_used_for_scaling = 50.;
        let temp_limit_upper = 310.;

        let result = heat_pump.energy_output_limited(
            energy_output_required,
            Some(temp_output),
            temp_used_for_scaling,
            temp_limit_upper,
        );

        assert_relative_eq!(result, 1.4444444444444444);

        let energy_output_required = 1.5;
        let temp_output = 320.;
        let temp_used_for_scaling = 320.;
        let temp_limit_upper = 310.;

        let result = heat_pump.energy_output_limited(
            energy_output_required,
            Some(temp_output),
            temp_used_for_scaling,
            temp_limit_upper,
        );

        assert_relative_eq!(result, 1.5);

        let energy_output_required = 1.5;
        let temp_output = 310.;
        let temp_used_for_scaling = 320.;
        let temp_limit_upper = 300.;

        let result = heat_pump.energy_output_limited(
            energy_output_required,
            Some(temp_output),
            temp_used_for_scaling,
            temp_limit_upper,
        );

        assert_relative_eq!(result, 0.);
    }

    #[rstest]
    /// Check if backup heater is available or still in delay period
    fn test_backup_heater_delay_time_elapsed(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
        let heat_pump = Arc::from(Mutex::from(heat_pump));
        let service_name = "service_backupheater";

        let _ = HeatPump::create_service_connection(heat_pump.clone(), service_name);

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            let _ = heat_pump.lock().demand_energy(
                service_name,
                &HeatingServiceType::DomesticHotWaterRegular,
                500.,
                Some(330.),
                330.,
                340.,
                design_flow_temp_op_cond_k(55.),
                1560.,
                true,
                t_it,
                Some(TempSpreadCorrectionArg::Float(1.)),
                None,
                None,
                None,
                None,
                None,
            );

            let result = heat_pump.lock().backup_heater_delay_time_elapsed();

            assert_eq!(result, [false, true][t_idx]);

            heat_pump.lock().timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    fn test_outside_operating_limits(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let result = heat_pump.outside_operating_limits(
            300.,
            simulation_time_for_heat_pump.iter().current_iteration(),
        );
        assert_eq!(result, false);

        let heat_pump_with_exhaust = create_heat_pump_with_exhaust(
            "aux_outside_operating_limit",
            Some(-1.),
            Some(10),
            external_conditions,
            simulation_time_for_heat_pump,
        );

        let result = heat_pump_with_exhaust.outside_operating_limits(
            300.,
            simulation_time_for_heat_pump.iter().current_iteration(),
        );
        assert_eq!(result, true)
    }

    #[rstest]
    fn test_outside_operating_limits_sinktypes(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );
        heat_pump.sink_type = HeatPumpSinkType::Water;

        assert!(!heat_pump.outside_operating_limits(
            300.,
            simulation_time_for_heat_pump.iter().current_iteration()
        ));
        assert!(heat_pump.outside_operating_limits(
            350.,
            simulation_time_for_heat_pump.iter().current_iteration()
        ));

        heat_pump.sink_type = HeatPumpSinkType::Air;

        assert!(!heat_pump.outside_operating_limits(
            350.,
            simulation_time_for_heat_pump.iter().current_iteration()
        ));
    }

    #[rstest]
    fn test_load_ratio_and_mode(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );
        // Test without modulating ctrl
        heat_pump.modulating_ctrl = false;

        let (load_ratio, load_ratio_continuous_min, hp_operating_in_onoff_mode) =
            heat_pump.load_ratio_and_mode(1., 1., 50.).unwrap();

        assert_relative_eq!(load_ratio, 1.);
        assert_relative_eq!(load_ratio_continuous_min, 1.);
        assert!(!hp_operating_in_onoff_mode);

        // Test with modulating ctrl
        heat_pump.modulating_ctrl = true;

        let (load_ratio, load_ratio_continuous_min, hp_operating_in_onoff_mode) =
            heat_pump.load_ratio_and_mode(0.6, 0.8, 50.).unwrap();

        assert_relative_eq!(load_ratio, 0.75);
        assert_relative_eq!(load_ratio_continuous_min, 0.35);
        assert!(!hp_operating_in_onoff_mode);

        // Test with zero time available
        let (load_ratio, load_ratio_continuous_min, hp_operating_in_onoff_mode) =
            heat_pump.load_ratio_and_mode(0., 0., 50.).unwrap();

        assert_relative_eq!(load_ratio, 0.);
        assert_relative_eq!(load_ratio_continuous_min, 0.35);
        assert!(!hp_operating_in_onoff_mode);

        // Test with zero time available and non-zero running time
        assert!(heat_pump.load_ratio_and_mode(0.2, 0., 50.).is_err());

        // Test HP with no test data for 55 degC design flow temp
        let heat_pump_no_55deg_data = create_heat_pump_with_energy_supply_heat_source(
            "arbitrary name for testing",
            external_conditions,
            simulation_time_for_heat_pump,
        );

        let (load_ratio, load_ratio_continuous_min, hp_operating_in_onoff_mode) =
            heat_pump_no_55deg_data
                .load_ratio_and_mode(0.2, 0.8, 45.)
                .unwrap();

        assert_relative_eq!(load_ratio, 0.25);
        assert_relative_eq!(load_ratio_continuous_min, 0.35);
        assert!(hp_operating_in_onoff_mode);
    }

    #[rstest]
    fn test_energy_input_compressor(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );
        let (energy_input_hp, compressor_power_min_load) = heat_pump.energy_input_compressor(
            true,
            false,
            false,
            8.4,
            8.4,
            Some(8.4),
            2.9,
            1.,
            1.,
            0.4,
            1560.,
            HeatingServiceType::DomesticHotWaterRegular,
        );

        assert_relative_eq!(energy_input_hp, 2.896551724137931);
        assert_relative_eq!(compressor_power_min_load, 1.1586206896551725);

        let (energy_input_hp, compressor_power_min_load) = heat_pump.energy_input_compressor(
            false,
            false,
            false,
            8.4,
            8.4,
            Some(8.4),
            2.9,
            1.,
            1.,
            0.4,
            1560.,
            HeatingServiceType::DomesticHotWaterRegular,
        );

        assert_relative_eq!(energy_input_hp, 0.);
        assert_relative_eq!(compressor_power_min_load, 1.1586206896551725);

        let (energy_input_hp, compressor_power_min_load) = heat_pump.energy_input_compressor(
            false,
            false,
            false,
            8.4,
            8.4,
            None,
            2.9,
            1.,
            1.,
            0.4,
            1560.,
            HeatingServiceType::DomesticHotWaterRegular,
        );

        assert_relative_eq!(energy_input_hp, 0.);
        assert_relative_eq!(compressor_power_min_load, 0.);

        let (energy_input_hp, compressor_power_min_load) = heat_pump.energy_input_compressor(
            true,
            false,
            true,
            8.4,
            8.4,
            Some(0.4),
            2.9,
            1.,
            1.,
            0.4,
            1560.,
            HeatingServiceType::Space,
        );

        assert_relative_eq!(energy_input_hp, 0.13793103448275862);
        assert_relative_eq!(compressor_power_min_load, 0.05517241379310345);

        let (energy_input_hp, compressor_power_min_load) = heat_pump.energy_input_compressor(
            true,
            false,
            true,
            8.4,
            8.4,
            Some(0.4),
            2.9,
            1.,
            1.,
            0.4,
            1560.,
            HeatingServiceType::DomesticHotWaterRegular,
        );

        assert_relative_eq!(energy_input_hp, 0.13793103448275862);
        assert_relative_eq!(compressor_power_min_load, 0.05517241379310345);

        heat_pump.sink_type = HeatPumpSinkType::Air;

        let (energy_input_hp, compressor_power_min_load) = heat_pump.energy_input_compressor(
            true,
            false,
            true,
            8.4,
            12.,
            Some(0.4),
            2.9,
            1.,
            1.,
            0.4,
            1560.,
            HeatingServiceType::Space,
        );

        assert_relative_eq!(energy_input_hp, 0.09655172413793105);
        assert_relative_eq!(compressor_power_min_load, 0.05517241379310345);

        let (energy_input_hp, compressor_power_min_load) = heat_pump.energy_input_compressor(
            true,
            false,
            true,
            8.4,
            8.4,
            Some(0.4),
            2.9,
            1.,
            1.,
            0.4,
            1560.,
            HeatingServiceType::DomesticHotWaterRegular,
        );

        assert_relative_eq!(energy_input_hp, 0.13793103448275862);
        assert_relative_eq!(compressor_power_min_load, 0.05517241379310345);
    }

    #[rstest]
    fn test_inadequate_capacity(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let energy_output_required = 5.0;
        let thermal_capacity_op_cond = 5.0;
        let temp_output = 310.0;
        let time_available = 1.0;
        let temp_return_feed = 315.0;
        let hybrid_boiler_service = None;
        let time_start = 0.;

        assert!(!heat_pump
            .inadequate_capacity(
                energy_output_required,
                thermal_capacity_op_cond,
                temp_output,
                time_available,
                time_start,
                temp_return_feed,
                hybrid_boiler_service,
                simulation_time_for_heat_pump.iter().current_iteration()
            )
            .unwrap());
    }

    #[rstest]
    fn test_inadequate_capacity_topup(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        heat_pump.set_canned_whether_delay_time_elapsed(Some(true));

        let energy_output_required = 5.0;
        let thermal_capacity_op_cond = 1.0;
        let temp_output = 343.;
        let time_available = 2.;
        let temp_return_feed = 313.;
        let hybrid_boiler_service = None;
        let time_start = 0.;

        let result = heat_pump
            .inadequate_capacity(
                energy_output_required,
                thermal_capacity_op_cond,
                temp_output,
                time_available,
                time_start,
                temp_return_feed,
                hybrid_boiler_service,
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        assert!(result);
    }

    #[rstest]
    fn test_inadequate_capacity_substitute(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            Some("Substitute"),
        );

        heat_pump.set_canned_whether_delay_time_elapsed(Some(true));

        let energy_output_required = 5.0;
        let thermal_capacity_op_cond = 1.0;
        let temp_output = 343.;
        let time_available = 2.;
        let temp_return_feed = 313.;
        let hybrid_boiler_service = None;
        let time_start = 0.;

        let result = heat_pump
            .inadequate_capacity(
                energy_output_required,
                thermal_capacity_op_cond,
                temp_output,
                time_available,
                time_start,
                temp_return_feed,
                hybrid_boiler_service,
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        assert!(result);
    }

    #[rstest]
    fn test_inadequate_capacity_no_delay_elapsed(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        heat_pump.set_canned_whether_delay_time_elapsed(Some(false));

        let energy_output_required = 5.0;
        let thermal_capacity_op_cond = 5.0;
        let temp_output = 343.;
        let time_available = 2.;
        let temp_return_feed = 313.;
        let hybrid_boiler_service = None;
        let time_start = 0.;

        let result = heat_pump
            .inadequate_capacity(
                energy_output_required,
                thermal_capacity_op_cond,
                temp_output,
                time_available,
                time_start,
                temp_return_feed,
                hybrid_boiler_service,
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        assert!(!result);
    }

    #[rstest]
    fn test_inadequate_capacity_insufficient_backup(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            Some("Substitute"),
        );

        heat_pump.set_canned_whether_delay_time_elapsed(Some(true));

        let energy_output_required = 5.0;
        let thermal_capacity_op_cond = 5.0;
        let temp_output = 343.;
        let time_available = 2.;
        let temp_return_feed = 313.;
        let hybrid_boiler_service = None;
        let time_start = 0.;

        let result = heat_pump
            .inadequate_capacity(
                energy_output_required,
                thermal_capacity_op_cond,
                temp_output,
                time_available,
                time_start,
                temp_return_feed,
                hybrid_boiler_service,
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        assert!(!result);
    }

    #[rstest]
    fn test_is_heat_pump_cost_effective_equal_cost(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let cop_op_cond = 3.;
        let boiler_eff = 0.9;

        let cost_schedule_hybrid_hp = json!({
           "cost_schedule_start_day": 0,
           "cost_schedule_time_series_step": 1,
           "cost_schedule_hp":
               {"main": [16, 20, 24, {"value": 16, "repeat": 5}]},
           "cost_schedule_boiler": {"main": [{"value": 4, "repeat": 8}]}
        });

        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: CostSchedule";

        let heat_pump_input = create_heat_pump_input_from_json(Some("Substitute"), None);

        let heat_pump_with_cost_schedule = create_heat_pump(
            heat_pump_input,
            energy_supply_conn_name_auxiliary,
            None, // energy_supply_heat_source: Option<Arc<RwLock<EnergySupply>>>,
            None, // boiler: Option<Arc<Mutex<Boiler>>>,
            None, // throughput_exhaust_air: Option<f64>,
            None, // temp_internal_air: Option<f64>,
            external_conditions,
            simulation_time_for_heat_pump,
            Some(cost_schedule_hybrid_hp),
        );

        let result = heat_pump_with_cost_schedule.is_heat_pump_cost_effective(
            cop_op_cond,
            boiler_eff,
            simulation_time_for_heat_pump.iter().current_iteration(),
        );

        assert!(!result);

        let cop_op_cond = 16.;
        let boiler_eff = 2.;

        let result = heat_pump_with_cost_schedule.is_heat_pump_cost_effective(
            cop_op_cond,
            boiler_eff,
            simulation_time_for_heat_pump.iter().current_iteration(),
        );

        assert!(result);
    }

    #[rstest]
    fn test_use_backup_heater_only(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let cop_op_cond = 3.0;
        let energy_output_required = 4.0;
        let thermal_capacity_op_cond = 4.0;
        let temp_output = 320.0;
        let time_available = 1.0;
        let temp_return_feed = 320.0;
        let time_start = 0.;

        let heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        assert!(!heat_pump
            .use_backup_heater_only(
                cop_op_cond,
                energy_output_required,
                thermal_capacity_op_cond,
                temp_output,
                time_available,
                time_start,
                temp_return_feed,
                None,
                None,
                simulation_time_for_heat_pump.iter().current_iteration()
            )
            .unwrap());

        let cost_schedule_hybrid_hp = json!({
           "cost_schedule_start_day": 0,
           "cost_schedule_time_series_step": 1,
           "cost_schedule_hp":
               {"main": [16, 20, 24, {"value": 16, "repeat": 5}]},
           "cost_schedule_boiler": {"main": [{"value": 4, "repeat": 8}]}
        });

        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: CostSchedule_1";

        let heat_pump_input = create_heat_pump_input_from_json(None, None);

        let heat_pump_with_cost_schedule = create_heat_pump(
            heat_pump_input,
            energy_supply_conn_name_auxiliary,
            None, // energy_supply_heat_source: Option<Arc<RwLock<EnergySupply>>>,
            None, // boiler: Option<Arc<Mutex<Boiler>>>,
            None, // throughput_exhaust_air: Option<f64>,
            None, // temp_internal_air: Option<f64>,
            external_conditions,
            simulation_time_for_heat_pump,
            Some(cost_schedule_hybrid_hp),
        );

        assert!(heat_pump_with_cost_schedule
            .use_backup_heater_only(
                cop_op_cond,
                energy_output_required,
                thermal_capacity_op_cond,
                temp_output,
                time_available,
                time_start,
                temp_return_feed,
                None,
                Some(3.0),
                simulation_time_for_heat_pump.iter().current_iteration()
            )
            .unwrap());
    }

    #[rstest]
    fn test_run_demand_energy_calc(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        // Test with hybrid_boiler_service and boiler_eff
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: boiler";

        let boiler = Arc::new(RwLock::new(create_boiler(
            external_conditions.clone(),
            energy_supply(simulation_time_for_heat_pump),
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));
        let control = SetpointTimeControl::new(
            vec![Some(21.), Some(22.)],
            0,
            1.0,
            Default::default(),
            Default::default(),
            1.0,
        );
        let boiler_service_space = Arc::new(Mutex::new(Boiler::create_service_space_heating(
            boiler.clone(),
            "service_boilerspace",
            Arc::new(Control::SetpointTime(control)),
        )));

        let heat_pump_input = create_heat_pump_input_from_json(None, None);

        let mut heat_pump_with_boiler = create_heat_pump(
            heat_pump_input,
            energy_supply_conn_name_auxiliary,
            None,
            Some(boiler),
            None,
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let expected_energy_calcs: [HeatPumpEnergyCalculation; 2] = [
            HeatPumpEnergyCalculation {
                service_name: "service_boilerspace".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(320.0),
                temp_source: 273.15,
                cop_op_cond: Some(3.395194725645021), // 3.39519472564502 in Python
                thermal_capacity_op_cond: Some(9.231749514150996),
                time_running_full_load: 0.10832182984028522,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 1.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 1.0,
                energy_heating_circ_pump: 0.0016248274476042782,
                energy_source_circ_pump: 0.0010832182984028523,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
            HeatPumpEnergyCalculation {
                service_name: "service_boilerspace".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(320.0),
                temp_source: 275.65,
                cop_op_cond: Some(3.515322632524056), // 3.5153226325240556 in Python
                thermal_capacity_op_cond: Some(9.487556781989081),
                time_running_full_load: 0.10540121371377431,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 1.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: Default::default(),
                energy_input_total: Default::default(),
                energy_delivered_total: 1.0,
                energy_heating_circ_pump: 0.0015810182057066145,
                energy_source_circ_pump: 0.0010540121371377432,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
        ];

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            assert_eq!(
                heat_pump_with_boiler
                    .run_demand_energy_calc(
                        "service_boilerspace",
                        &HeatingServiceType::DomesticHotWaterRegular,
                        1.,
                        Some(320.),
                        310.,
                        340.,
                        design_flow_temp_op_cond_k(55.),
                        1560.,
                        true,
                        t_it,
                        Some(TempSpreadCorrectionArg::Float(1.0)),
                        None,
                        Some(HybridBoilerService::Space(boiler_service_space.clone())),
                        Some(1.0),
                        Some(0.0),
                        None,
                        None,
                    )
                    .unwrap(),
                expected_energy_calcs[t_idx]
            );
        }

        // Check without boiler and service_on True
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let expected_energy_calcs: [HeatPumpEnergyCalculation; 2] = [
            HeatPumpEnergyCalculation {
                service_name: "service_no_boiler".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 273.15,
                cop_op_cond: Some(3.182442675928905), // 3.1824426759289044 in Python
                thermal_capacity_op_cond: Some(8.417674488123662),
                time_running_full_load: 0.11879765621857688,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 0.9999999999999999,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 0.9999999999999999,
                energy_heating_circ_pump: 0.001781964843278653,
                energy_source_circ_pump: 0.0011879765621857687,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
            HeatPumpEnergyCalculation {
                service_name: "service_no_boiler".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 275.65,
                cop_op_cond: Some(3.1971341044162673), // 3.197134104416267 in Python
                thermal_capacity_op_cond: Some(8.650924134797519),
                time_running_full_load: 0.11559458670751663,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 1.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 1.0,
                energy_heating_circ_pump: 0.0017339188006127494,
                energy_source_circ_pump: 0.0011559458670751662,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
        ];

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            assert_eq!(
                heat_pump
                    .run_demand_energy_calc(
                        "service_no_boiler",
                        &HeatingServiceType::DomesticHotWaterRegular,
                        1.,
                        Some(330.),
                        330.,
                        340.,
                        design_flow_temp_op_cond_k(55.),
                        1560.,
                        true,
                        t_it,
                        Some(TempSpreadCorrectionArg::Float(1.0)),
                        None,
                        None,
                        None,
                        Some(0.0),
                        None,
                        None,
                    )
                    .unwrap(),
                expected_energy_calcs[t_idx]
            );
        }

        let heat_modcontrol_input = serde_json::from_value(json!(
            {
                "type": "HeatPump",
                "EnergySupply": "mains_gas",
                "source_type": "OutsideAir",
                "sink_type": "Water",
                "backup_ctrl_type": "Substitute",
                "time_delay_backup": 2.0,
                "modulating_control": false,
                "min_modulation_rate_35": 0.35,
                "min_modulation_rate_55": 0.4,
                "time_constant_onoff_operation": 140,
                "temp_return_feed_max": 70.0,
                "temp_lower_operating_limit": -5.0,
                "min_temp_diff_flow_return_for_hp_to_operate": 0.0,
                "var_flow_temp_ctrl_during_test": true,
                "power_heating_circ_pump": 0.015,
                "power_source_circ_pump": 0.01,
                "power_standby": 0.015,
                "power_crankcase_heater": 0.01,
                "power_off": 0.015,
                "power_max_backup": 3.0,
                "test_data_EN14825": [{"test_letter": "A",
                               "capacity": 8.4,
                               "cop": 4.6,
                               "design_flow_temp": 35,
                               "temp_outlet": 34,
                               "temp_source": 0,
                               "temp_test": -7},
                               {"test_letter": "B",
                                 "capacity": 8.3,
                                 "cop": 4.9,
                                 "design_flow_temp": 35,
                                 "temp_outlet": 30,
                                 "temp_source": 0,
                                 "temp_test": 2},
                               {"test_letter": "C",
                                "capacity": 8.3,
                                "cop": 5.1,
                                "design_flow_temp": 35,
                                "temp_outlet": 27,
                                "temp_source": 0,
                                "temp_test": 7},
                               {"test_letter": "D",
                                "capacity": 8.2,
                                "cop": 5.4,
                                "design_flow_temp": 35,
                                "temp_outlet": 24,
                                "temp_source": 0,
                                "temp_test": 12},
                                {"test_letter": "F",
                                 "capacity": 8.4,
                                 "cop": 4.6,
                                 "design_flow_temp": 35,
                                 "temp_outlet": 34,
                                 "temp_source": 0,
                                 "temp_test": -7}
                               ]
                }
        ))
        .unwrap();

        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: modulating_control";

        let mut heat_pump_mod_ctrl = create_heat_pump(
            heat_modcontrol_input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            None,
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let expected_modcontrol_calcs = [
            HeatPumpEnergyCalculation {
                service_name: "service_backup_modctrl".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 273.15,
                cop_op_cond: Some(3.9929681311183423),
                thermal_capacity_op_cond: Some(8.857000000000003),
                time_running_full_load: 0.11290504685559441,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 1.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 1.0,
                energy_heating_circ_pump: 0.0016935757028339162,
                energy_source_circ_pump: 0.0011290504685559442,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
            HeatPumpEnergyCalculation {
                service_name: "service_backup_modctrl".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 275.65,
                cop_op_cond: Some(4.089173762601537),
                thermal_capacity_op_cond: Some(8.807000000000002),
                time_running_full_load: 0.1135460429204042,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 1.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 1.0,
                energy_heating_circ_pump: 0.001703190643806063,
                energy_source_circ_pump: 0.001135460429204042,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
        ];

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            assert_eq!(
                heat_pump_mod_ctrl
                    .run_demand_energy_calc(
                        "service_backup_modctrl",
                        &HeatingServiceType::DomesticHotWaterRegular,
                        1.0,
                        Some(330.0),
                        330.0,
                        340.0,
                        design_flow_temp_op_cond_k(55.),
                        1560.,
                        true,
                        t_it,
                        Some(TempSpreadCorrectionArg::Float(1.0)),
                        None,
                        None,
                        None,
                        Some(0.0),
                        None,
                        None,
                    )
                    .unwrap(),
                expected_modcontrol_calcs[t_idx]
            );
        }

        // Check the results with service off and more energy_output_required
        let expected_energy_calcs: [HeatPumpEnergyCalculation; 2] = [
            HeatPumpEnergyCalculation {
                service_name: "service_energy_output_required".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: false,
                energy_output_required: 50.,
                temp_output: Some(330.0),
                temp_source: 273.15,
                cop_op_cond: Some(3.182442675928905), // 3.1824426759289044 in Python
                thermal_capacity_op_cond: Some(8.417674488123662),
                time_running_full_load: 0.0,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 0.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 0.0,
                energy_heating_circ_pump: 0.0,
                energy_source_circ_pump: 0.0,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
            HeatPumpEnergyCalculation {
                service_name: "service_energy_output_required".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: false,
                energy_output_required: 50.,
                temp_output: Some(330.0),
                temp_source: 275.65,
                cop_op_cond: Some(3.1971341044162673), // 3.197134104416267 in Python
                thermal_capacity_op_cond: Some(8.650924134797519),
                time_running_full_load: 0.0,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 0.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 0.0,
                energy_heating_circ_pump: 0.0,
                energy_source_circ_pump: 0.0,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
        ];

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            assert_eq!(
                heat_pump
                    .run_demand_energy_calc(
                        "service_energy_output_required",
                        &HeatingServiceType::DomesticHotWaterRegular,
                        50.,
                        Some(330.),
                        330.,
                        340.,
                        design_flow_temp_op_cond_k(55.),
                        1560.,
                        false,
                        t_it,
                        Some(TempSpreadCorrectionArg::Float(1.)),
                        None,
                        None,
                        None,
                        Some(0.0),
                        None,
                        None,
                    )
                    .unwrap(),
                expected_energy_calcs[t_idx]
            );
        }

        // Check service off with boiler

        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: boiler_and_sevice_off";

        let boiler = Arc::new(RwLock::new(create_boiler(
            external_conditions.clone(),
            energy_supply(simulation_time_for_heat_pump),
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));

        let ctrl = Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(21.0), Some(22.0)],
            0,
            1.0,
            Default::default(),
            Default::default(),
            1.0,
        ));

        let boiler_service_space = Arc::new(Mutex::new(Boiler::create_service_space_heating(
            boiler.clone(),
            "service_boilerspace_service_off",
            ctrl.into(),
        )));

        let heat_pump_input = create_heat_pump_input_from_json(None, None);

        let mut heat_pump_with_boiler = create_heat_pump(
            heat_pump_input,
            energy_supply_conn_name_auxiliary,
            None,
            Some(boiler),
            None,
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        let expected_off_with_boiler_calcs = [
            HeatPumpEnergyCalculation {
                service_name: "service_boilerspace_service_off".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: false,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 273.15,
                cop_op_cond: Some(3.182442675928905), // 3.1824426759289044 in Python
                thermal_capacity_op_cond: Some(8.417674488123662),
                time_running_full_load: 0.0,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 0.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 0.0,
                energy_heating_circ_pump: 0.0,
                energy_source_circ_pump: 0.0,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
            HeatPumpEnergyCalculation {
                service_name: "service_boilerspace_service_off".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: false,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 275.65,
                cop_op_cond: Some(3.1971341044162673), // 3.197134104416267 in Python
                thermal_capacity_op_cond: Some(8.650924134797519),
                time_running_full_load: 0.0,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 0.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 0.0,
                energy_heating_circ_pump: 0.0,
                energy_source_circ_pump: 0.0,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
        ];

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            assert_eq!(
                heat_pump_with_boiler
                    .run_demand_energy_calc(
                        "service_boilerspace_service_off",
                        &HeatingServiceType::DomesticHotWaterRegular,
                        1.0,
                        Some(330.0),
                        330.0,
                        340.0,
                        design_flow_temp_op_cond_k(55.),
                        1560.,
                        false,
                        t_it,
                        Some(TempSpreadCorrectionArg::Float(1.0)),
                        None,
                        Some(HybridBoilerService::Space(boiler_service_space.clone())),
                        Some(1.0),
                        Some(0.0),
                        None,
                        None,
                    )
                    .unwrap(),
                expected_off_with_boiler_calcs[t_idx]
            );
        }

        // Check with backup_ctrl_type Substitute
        let heat_pump_with_sub_input = serde_json::from_value(json!(
            {"type": "HeatPump",
                                    "EnergySupply": "mains_gas",
                                    "source_type": "OutsideAir",
                                    "sink_type": "Water",
                                    "backup_ctrl_type": "Substitute",
                                    "time_delay_backup": 2.0,
                                    "modulating_control": true,
                                    "min_modulation_rate_35": 0.35,
                                    "min_modulation_rate_55": 0.4,
                                    "time_constant_onoff_operation": 140,
                                    "temp_return_feed_max": 70.0,
                                    "temp_lower_operating_limit": -5.0,
                                    "min_temp_diff_flow_return_for_hp_to_operate": 0.0,
                                    "var_flow_temp_ctrl_during_test": true,
                                    "power_heating_circ_pump": 0.015,
                                    "power_source_circ_pump": 0.01,
                                    "power_standby": 0.015,
                                    "power_crankcase_heater": 0.01,
                                    "power_off": 0.015,
                                    "power_max_backup": 3.0,
                                    "test_data_EN14825": [{"test_letter": "A",
                                                   "capacity": 8.4,
                                                   "cop": 4.6,
                                                   "design_flow_temp": 35,
                                                   "temp_outlet": 34,
                                                   "temp_source": 0,
                                                   "temp_test": -7},
                                                   {"test_letter": "B",
                                                     "capacity": 8.3,
                                                     "cop": 4.9,
                                                     "design_flow_temp": 35,
                                                     "temp_outlet": 30,
                                                     "temp_source": 0,
                                                     "temp_test": 2},
                                                   {"test_letter": "C",
                                                    "capacity": 8.3,
                                                    "cop": 5.1,
                                                    "design_flow_temp": 35,
                                                    "temp_outlet": 27,
                                                    "temp_source": 0,
                                                    "temp_test": 7},
                                                   {"test_letter": "D",
                                                    "capacity": 8.2,
                                                    "cop": 5.4,
                                                    "design_flow_temp": 35,
                                                    "temp_outlet": 24,
                                                    "temp_source": 0,
                                                    "temp_test": 12},
                                                    {"test_letter": "F",
                                                     "capacity": 8.4,
                                                     "cop": 4.6,
                                                     "design_flow_temp": 35,
                                                     "temp_outlet": 34,
                                                     "temp_source": 0,
                                                     "temp_test": -7}
                                                   ]
                                    }
        ))
        .unwrap();

        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: backup";

        let mut heat_pump_backup = create_heat_pump(
            heat_pump_with_sub_input,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        let expected_energy_calcs: [HeatPumpEnergyCalculation; 2] = [
            HeatPumpEnergyCalculation {
                service_name: "service_backup_substitute".into(),
                service_type: HeatingServiceType::Space,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 273.15,
                cop_op_cond: Some(3.9929681311183423),
                thermal_capacity_op_cond: Some(6.773123981338176),
                time_running_full_load: 0.1476423586450323,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 1.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 1.0,
                energy_heating_circ_pump: 0.0022146353796754846,
                energy_source_circ_pump: 0.0014764235864503231,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
            HeatPumpEnergyCalculation {
                service_name: "service_backup_substitute".into(),
                service_type: HeatingServiceType::Space,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 275.65,
                cop_op_cond: Some(4.089173762601537),
                thermal_capacity_op_cond: Some(6.9608039370972525),
                time_running_full_load: 0.1436615668300253,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 1.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 1.0,
                energy_heating_circ_pump: 0.002154923502450379,
                energy_source_circ_pump: 0.0014366156683002528,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            },
        ];

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            assert_eq!(
                heat_pump_backup
                    .run_demand_energy_calc(
                        "service_backup_substitute",
                        &HeatingServiceType::Space,
                        1.0,
                        Some(330.0),
                        330.0,
                        340.0,
                        design_flow_temp_op_cond_k(55.),
                        1560.,
                        true,
                        t_it,
                        Some(TempSpreadCorrectionArg::Float(1.0)),
                        None,
                        None,
                        None,
                        Some(0.0),
                        None,
                        None,
                    )
                    .unwrap(),
                expected_energy_calcs[t_idx]
            )
        }
    }

    #[rstest]
    fn test_run_demand_energy_calc_other(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
        mut buffer_tank: BufferTank,
    ) {
        let emitters_data_for_buffer_tank = BufferTankEmittersDataWithResult {
            data: {
                BufferTankEmittersData {
                    temp_emitter_req: 43.,
                    power_req_from_buffer_tank: 6.,
                    design_flow_temp: 55.,
                    target_flow_temp: 48.,
                    temp_rm_prev: 22.,
                    variable_flow: true,
                    temp_diff_emit_design: 10.,
                    min_flow_rate: 0.05,
                    max_flow_rate: 0.3,
                }
            },
            result: {
                BufferTankServiceResult {
                    _service_name: Default::default(),
                    _power_req_from_buffer_tank: Default::default(),
                    _temp_emitter_req: Default::default(),
                    _buffer_emitter_circ_flow_rate: Default::default(),
                    flow_temp_increase_due_to_buffer: 10.,
                    pump_power_at_flow_rate: 20.,
                    heat_loss_buffer_kwh: 30.,
                }
            },
        };
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hp";
        let boiler = Arc::new(RwLock::new(create_boiler(
            external_conditions.clone(),
            energy_supply(simulation_time_for_heat_pump),
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));
        let control = create_setpoint_time_control(vec![Some(21.), Some(22.)]);

        let boiler_service_space = Arc::new(Mutex::new(Boiler::create_service_space_heating(
            boiler.clone(),
            "service_boilerspace",
            Arc::new(control),
        )));
        let input = create_heat_pump_input_from_json(None, None);
        let mut heat_pump_with_boiler = create_heat_pump(
            input,
            energy_supply_conn_name_auxiliary,
            None,
            Some(boiler),
            None,
            None,
            external_conditions.clone(),
            simulation_time_for_heat_pump,
            None,
        );

        // Test with ServiceType Water
        let results = heat_pump_with_boiler
            .run_demand_energy_calc(
                "service_boilerspace",
                &HeatingServiceType::DomesticHotWaterRegular,
                1.,
                Some(320.0),
                310.0,
                340.0,
                design_flow_temp_op_cond_k(55.),
                1560.,
                true,
                simulation_time_for_heat_pump.iter().current_iteration(),
                Some(TempSpreadCorrectionArg::Float(1.0)),
                None,
                Some(HybridBoilerService::Space(boiler_service_space.clone())),
                Some(1.),
                Some(0.0),
                None,
                None,
            )
            .unwrap();

        let expected_results = HeatPumpEnergyCalculation {
            service_name: "service_boilerspace".into(),
            service_type: HeatingServiceType::DomesticHotWaterRegular,
            service_on: true,
            energy_output_required: 1.0,
            temp_output: Some(320.),
            temp_source: 273.15,
            cop_op_cond: Some(3.395194725645021), // 3.39519472564502 in Python
            thermal_capacity_op_cond: Some(9.231749514150996),
            time_running_full_load: 0.10832182984028522,
            time_running_part_load: None,
            time_available_for_capacity_calc: 1.0,
            time_available_for_load_ratio_calc: 1.0,
            time_constant_for_service: 1560.,
            compressor_power_min_load: Default::default(),
            load_ratio_continuous_min: Default::default(),
            load_ratio: Default::default(),
            use_backup_heater_only: false,
            hp_operating_in_onoff_mode: Default::default(),
            energy_input_hp: Default::default(),
            energy_delivered_hp: 1.0,
            energy_input_backup: 0.0,
            energy_delivered_backup: 0.0,
            energy_input_total: Default::default(),
            energy_delivered_total: 1.0,
            energy_heating_circ_pump: 0.0016248274476042782,
            energy_source_circ_pump: 0.0010832182984028523,
            energy_output_required_boiler: 0.0,
            energy_heating_warm_air_fan: 0.0,
            energy_output_delivered_boiler: None,
        };

        assert_eq!(results, expected_results);

        // Test with service_on set to false
        let results = heat_pump_with_boiler
            .run_demand_energy_calc(
                "service_boilerspace",
                &HeatingServiceType::DomesticHotWaterRegular,
                1.,
                None,
                310.0,
                340.0,
                design_flow_temp_op_cond_k(55.),
                1560.,
                false,
                simulation_time_for_heat_pump.iter().current_iteration(),
                Some(TempSpreadCorrectionArg::Float(1.0)),
                None,
                Some(HybridBoilerService::Space(boiler_service_space.clone())),
                Some(1.),
                Some(0.0),
                None,
                None,
            )
            .unwrap();

        let expected_results = HeatPumpEnergyCalculation {
            service_name: "service_boilerspace".into(),
            service_type: HeatingServiceType::DomesticHotWaterRegular,
            service_on: false,
            energy_output_required: 1.0,
            temp_output: None,
            temp_source: 273.15,
            cop_op_cond: None,
            thermal_capacity_op_cond: None,
            time_running_full_load: 0.,
            time_running_part_load: None,
            time_available_for_capacity_calc: 1.0,
            time_available_for_load_ratio_calc: 1.0,
            time_constant_for_service: 1560.,
            compressor_power_min_load: Default::default(),
            load_ratio_continuous_min: Default::default(),
            load_ratio: Default::default(),
            use_backup_heater_only: false,
            hp_operating_in_onoff_mode: Default::default(),
            energy_input_hp: Default::default(),
            energy_delivered_hp: 0.0,
            energy_input_backup: 0.0,
            energy_delivered_backup: 0.0,
            energy_input_total: Default::default(),
            energy_delivered_total: 0.0,
            energy_heating_circ_pump: 0.,
            energy_source_circ_pump: 0.,
            energy_output_required_boiler: 0.0,
            energy_heating_warm_air_fan: 0.0,
            energy_output_delivered_boiler: None,
        };

        assert_eq!(results, expected_results);

        // Test with HeatPumpBackupControlType.SUBSTITUTE and HeatPumpSinkType.AIR
        heat_pump_with_boiler.backup_ctrl = HeatPumpBackupControlType::Substitute;
        heat_pump_with_boiler.sink_type = HeatPumpSinkType::Air;

        let results = heat_pump_with_boiler
            .run_demand_energy_calc(
                "service_boilerspace",
                &HeatingServiceType::Space,
                1.,
                Some(320.),
                310.0,
                340.0,
                design_flow_temp_op_cond_k(55.),
                1560.,
                true,
                simulation_time_for_heat_pump.iter().current_iteration(),
                Some(TempSpreadCorrectionArg::Float(1.0)),
                None,
                Some(HybridBoilerService::Space(boiler_service_space.clone())),
                Some(1.),
                Some(0.0),
                None,
                None,
            )
            .unwrap();

        let expected_results = HeatPumpEnergyCalculation {
            service_name: "service_boilerspace".into(),
            service_type: HeatingServiceType::Space,
            service_on: true,
            energy_output_required: 1.0,
            temp_output: Some(320.),
            temp_source: 273.15,
            cop_op_cond: Some(3.395194725645021), // 3.39519472564502 in Python
            thermal_capacity_op_cond: Some(9.231749514150996),
            time_running_full_load: 0.10832182984028522,
            time_running_part_load: None,
            time_available_for_capacity_calc: 1.0,
            time_available_for_load_ratio_calc: 1.0,
            time_constant_for_service: 1560.,
            compressor_power_min_load: Default::default(),
            load_ratio_continuous_min: Default::default(),
            load_ratio: Default::default(),
            use_backup_heater_only: false,
            hp_operating_in_onoff_mode: Default::default(),
            energy_input_hp: Default::default(),
            energy_delivered_hp: 1.0,
            energy_input_backup: 0.0,
            energy_delivered_backup: 0.0,
            energy_input_total: Default::default(),
            energy_delivered_total: 1.0,
            energy_heating_circ_pump: 0.,
            energy_source_circ_pump: 0.0010832182984028523,
            energy_output_required_boiler: 0.0,
            energy_heating_warm_air_fan: 0.0,
            energy_output_delivered_boiler: None,
        };

        assert_eq!(results, expected_results);

        buffer_tank.set_track_buffer_loss(5.);
        // Python sets buffer_loss return value to 5 (buffer_tank.get_buffer_loss.return_value = 5)
        assert_eq!(buffer_tank.buffer_loss(), 5.);

        heat_pump_with_boiler.buffer_tank = Some(buffer_tank.clone());

        // Test with a buffer tank and emitters_data_for_buffer_tank
        heat_pump_with_boiler.time_running_continuous = 1.;

        let results = heat_pump_with_boiler
            .run_demand_energy_calc(
                "service_boilerspace",
                &HeatingServiceType::Space,
                1.,
                Some(320.),
                310.0,
                340.0,
                design_flow_temp_op_cond_k(55.),
                1560.,
                true,
                simulation_time_for_heat_pump.iter().current_iteration(),
                Some(TempSpreadCorrectionArg::Float(1.0)),
                None,
                Some(HybridBoilerService::Space(boiler_service_space.clone())),
                Some(1.),
                Some(0.0),
                None,
                Some(emitters_data_for_buffer_tank.clone()),
            )
            .unwrap();

        let expected_results = HeatPumpEnergyCalculation {
            service_name: "service_boilerspace".into(),
            service_type: HeatingServiceType::Space,
            service_on: true,
            energy_output_required: 36.0,
            temp_output: Some(330.),
            temp_source: 273.15,
            cop_op_cond: Some(3.182442675928905), // 3.1824426759289044 in Python
            thermal_capacity_op_cond: Some(8.417674488123662),
            time_running_full_load: 0.,
            time_running_part_load: None,
            time_available_for_capacity_calc: 1.0,
            time_available_for_load_ratio_calc: 1.0,
            time_constant_for_service: 1560.,
            compressor_power_min_load: Default::default(),
            load_ratio_continuous_min: Default::default(),
            load_ratio: Default::default(),
            use_backup_heater_only: true,
            hp_operating_in_onoff_mode: Default::default(),
            energy_input_hp: Default::default(),
            energy_delivered_hp: 0.0,
            energy_input_backup: 0.0,
            energy_delivered_backup: 0.0,
            energy_input_total: Default::default(),
            energy_delivered_total: 0.0,
            energy_heating_circ_pump: 0.,
            energy_source_circ_pump: 0.,
            energy_output_required_boiler: 24.0,
            energy_heating_warm_air_fan: 0.0,
            energy_output_delivered_boiler: None,
        };

        assert_eq!(results, expected_results);

        // Test with a buffer tank and emitters_data_for_buffer_tank, and service_on set to false
        let results = heat_pump_with_boiler
            .run_demand_energy_calc(
                "service_boilerspace",
                &HeatingServiceType::Space,
                1.,
                Some(320.),
                310.0,
                340.0,
                design_flow_temp_op_cond_k(55.),
                1560.,
                false,
                simulation_time_for_heat_pump.iter().current_iteration(),
                Some(TempSpreadCorrectionArg::Float(1.0)),
                None,
                Some(HybridBoilerService::Space(boiler_service_space.clone())),
                Some(1.),
                Some(0.0),
                None,
                Some(emitters_data_for_buffer_tank.clone()),
            )
            .unwrap();

        let expected_results = HeatPumpEnergyCalculation {
            service_name: "service_boilerspace".into(),
            service_type: HeatingServiceType::Space,
            service_on: false,
            energy_output_required: 1.0,
            temp_output: Some(330.),
            temp_source: 273.15,
            cop_op_cond: Some(3.182442675928905), // 3.1824426759289044 in Python
            thermal_capacity_op_cond: Some(8.417674488123662),
            time_running_full_load: 0.,
            time_running_part_load: None,
            time_available_for_capacity_calc: 1.0,
            time_available_for_load_ratio_calc: 1.0,
            time_constant_for_service: 1560.,
            compressor_power_min_load: Default::default(),
            load_ratio_continuous_min: Default::default(),
            load_ratio: Default::default(),
            use_backup_heater_only: false,
            hp_operating_in_onoff_mode: Default::default(),
            energy_input_hp: Default::default(),
            energy_delivered_hp: 0.0,
            energy_input_backup: 0.0,
            energy_delivered_backup: 0.0,
            energy_input_total: Default::default(),
            energy_delivered_total: 0.0,
            energy_heating_circ_pump: 0.,
            energy_source_circ_pump: 0.,
            energy_output_required_boiler: 0.0,
            energy_heating_warm_air_fan: 0.0,
            energy_output_delivered_boiler: None,
        };

        assert_eq!(results, expected_results);

        // Test with a buffer tank and energy delivered greater than buffer heat loss
        heat_pump_with_boiler.buffer_tank = Some(buffer_tank.clone());

        let results = heat_pump_with_boiler
            .run_demand_energy_calc(
                "TEST",
                &HeatingServiceType::Space,
                1.,
                Some(100.),
                310.0,
                340.0,
                design_flow_temp_op_cond_k(55.),
                1560.,
                true,
                simulation_time_for_heat_pump.iter().current_iteration(),
                Some(TempSpreadCorrectionArg::Float(1.0)),
                None,
                Some(HybridBoilerService::Space(boiler_service_space.clone())),
                Some(1.),
                Some(0.0),
                None,
                Some(emitters_data_for_buffer_tank),
            )
            .unwrap();

        let expected_results = HeatPumpEnergyCalculation {
            service_name: "TEST".into(),
            service_type: HeatingServiceType::Space,
            service_on: true,
            energy_output_required: 36.0,
            temp_output: Some(110.),
            temp_source: 273.15,
            cop_op_cond: Some(1.),
            thermal_capacity_op_cond: Some(227.27721117933888),
            time_running_full_load: 0.1583968749581025,
            time_running_part_load: None,
            time_available_for_capacity_calc: 1.0,
            time_available_for_load_ratio_calc: 1.0,
            time_constant_for_service: 1560.,
            compressor_power_min_load: Default::default(),
            load_ratio_continuous_min: Default::default(),
            load_ratio: Default::default(),
            use_backup_heater_only: false,
            hp_operating_in_onoff_mode: Default::default(),
            energy_input_hp: Default::default(),
            energy_delivered_hp: 36.0,
            energy_input_backup: 0.0,
            energy_delivered_backup: 0.0,
            energy_input_total: Default::default(),
            energy_delivered_total: 1.0,
            energy_heating_circ_pump: 0.,
            energy_source_circ_pump: 0.0015839687495810251,
            energy_output_required_boiler: 0.0,
            energy_heating_warm_air_fan: 0.0,
            energy_output_delivered_boiler: None,
        };

        assert_eq!(results, expected_results);
    }

    #[rstest]
    fn test_demand_energy(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
        energy_supply: EnergySupply,
    ) {
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: boiler";

        let boiler = Arc::new(RwLock::new(create_boiler(
            external_conditions.clone(),
            energy_supply,
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));

        let ctrl = Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(21.0), Some(22.0)],
            0,
            1.0,
            Default::default(),
            Default::default(),
            1.0,
        ));

        let boiler_service_space = Arc::new(Mutex::new(Boiler::create_service_space_heating(
            boiler.clone(),
            "service_boilerspace",
            Arc::new(ctrl),
        )));

        let heat_pump_input = create_heat_pump_input_from_json(None, None);

        let heat_pump_with_boiler = Arc::new(Mutex::new(create_heat_pump(
            heat_pump_input,
            energy_supply_conn_name_auxiliary,
            None,
            Some(boiler),
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        )));

        HeatPump::create_service_connection(
            heat_pump_with_boiler.clone(),
            "service_boiler_demand_energy",
        )
        .unwrap();

        let mut heat_pump_with_boiler = heat_pump_with_boiler.lock();

        for t_it in simulation_time_for_heat_pump.iter() {
            assert_relative_eq!(
                heat_pump_with_boiler
                    .demand_energy(
                        "service_boiler_demand_energy",
                        &HeatingServiceType::DomesticHotWaterRegular,
                        1.0,
                        Some(330.0),
                        330.0,
                        340.0,
                        design_flow_temp_op_cond_k(55.),
                        1560.,
                        true,
                        t_it,
                        Some(TempSpreadCorrectionArg::Float(1.0)),
                        None,
                        None,
                        Some(HybridBoilerService::Space(boiler_service_space.clone())),
                        None,
                        None,
                    )
                    .unwrap(),
                1.0
            );
        }

        assert_eq!(
            *heat_pump_with_boiler.service_results.read().deref(),
            vec![
                ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
                    service_name: "service_boiler_demand_energy".into(),
                    service_type: HeatingServiceType::DomesticHotWaterRegular,
                    service_on: true,
                    energy_output_required: 1.0,
                    temp_output: Some(330.0),
                    temp_source: 273.15,
                    cop_op_cond: Some(3.182442675928905), // 3.1824426759289044 in Python
                    thermal_capacity_op_cond: Some(8.417674488123662),
                    time_running_full_load: 0.11879765621857688,
                    time_running_part_load: None,
                    time_available_for_capacity_calc: 1.,
                    time_available_for_load_ratio_calc: 1.,
                    time_constant_for_service: 1560.,
                    compressor_power_min_load: Default::default(),
                    load_ratio_continuous_min: Default::default(),
                    load_ratio: Default::default(),
                    use_backup_heater_only: false,
                    hp_operating_in_onoff_mode: Default::default(),
                    energy_input_hp: Default::default(),
                    energy_delivered_hp: 0.9999999999999999,
                    energy_input_backup: 0.0,
                    energy_delivered_backup: 0.0,
                    energy_input_total: Default::default(),
                    energy_delivered_total: 0.9999999999999999,
                    energy_heating_circ_pump: 0.001781964843278653,
                    energy_source_circ_pump: 0.0011879765621857687,
                    energy_output_required_boiler: 0.0,
                    energy_heating_warm_air_fan: 0.,
                    energy_output_delivered_boiler: Some(0.0)
                })),
                ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
                    service_name: "service_boiler_demand_energy".into(),
                    service_type: HeatingServiceType::DomesticHotWaterRegular,
                    service_on: true,
                    energy_output_required: 1.0,
                    temp_output: Some(330.0),
                    temp_source: 275.65,
                    cop_op_cond: Some(3.1971341044162673), // 3.197134104416267 in Python
                    thermal_capacity_op_cond: Some(8.650924134797519),
                    time_running_full_load: 0.11559458670751663,
                    time_running_part_load: None,
                    time_available_for_capacity_calc: 0.8812023437814231,
                    time_available_for_load_ratio_calc: 0.8812023437814231,
                    time_constant_for_service: 1560.,
                    compressor_power_min_load: Default::default(),
                    load_ratio_continuous_min: Default::default(),
                    load_ratio: Default::default(),
                    use_backup_heater_only: false,
                    hp_operating_in_onoff_mode: Default::default(),
                    energy_input_hp: Default::default(),
                    energy_delivered_hp: 1.0,
                    energy_input_backup: 0.0,
                    energy_delivered_backup: 0.0,
                    energy_input_total: Default::default(),
                    energy_delivered_total: 1.0,
                    energy_heating_circ_pump: 0.0017339188006127494,
                    energy_source_circ_pump: 0.0011559458670751662,
                    energy_output_required_boiler: 0.0,
                    energy_heating_warm_air_fan: 0.,
                    energy_output_delivered_boiler: Some(0.0)
                })),
            ]
        );
    }

    #[rstest]
    fn test_demand_energy_hybrid_boiler_service(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
        energy_supply: EnergySupply,
    ) {
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: boiler";

        let boiler = Arc::new(RwLock::new(create_boiler(
            external_conditions.clone(),
            energy_supply,
            simulation_time_for_heat_pump,
            energy_supply_conn_name_auxiliary,
        )));

        let ctrl = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(21.0), Some(22.0)],
            0,
            1.0,
            Default::default(),
            Default::default(),
            1.0,
        )));

        let boiler_service_water = Arc::new(Mutex::new(
            Boiler::create_service_hot_water_regular(
                boiler.clone(),
                "service_water",
                ctrl.clone(),
                ctrl,
            )
            .unwrap(),
        ));

        let heat_pump_input = create_heat_pump_input_from_json(Some("Substitute"), None);

        let heat_pump_with_boiler = Arc::new(Mutex::new(create_heat_pump(
            heat_pump_input,
            energy_supply_conn_name_auxiliary,
            None,
            Some(boiler),
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        )));

        HeatPump::create_service_connection(
            heat_pump_with_boiler.clone(),
            "service_boiler_demand_energy",
        )
        .unwrap();

        let results = heat_pump_with_boiler
            .lock()
            .demand_energy(
                "service_boiler_demand_energy",
                &HeatingServiceType::DomesticHotWaterRegular,
                1.0,
                Some(330.0),
                330.0,
                340.0,
                design_flow_temp_op_cond_k(55.),
                1560.,
                true,
                simulation_time_for_heat_pump.iter().current_iteration(),
                Some(TempSpreadCorrectionArg::Float(1.0)),
                None,
                None,
                Some(HybridBoilerService::Regular(boiler_service_water.clone())),
                None,
                None,
            )
            .unwrap();

        assert_relative_eq!(results, 1.);
        assert_eq!(
            heat_pump_with_boiler
                .lock()
                .total_time_running_current_timestep_full_load,
            0.11879765621857688
        );
    }

    #[rstest]
    fn test_throughput_factor(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
        heat_pump.total_time_running_current_timestep_full_load = 2.;
        heat_pump.overvent_ratio = 0.8;

        assert_relative_eq!(heat_pump.throughput_factor(), 0.6);

        heat_pump.total_time_running_current_timestep_full_load = 3.;
        heat_pump.overvent_ratio = 0.7;

        assert_relative_eq!(heat_pump.throughput_factor(), 0.09999999999999964);
    }

    #[rstest]
    fn test_running_time_throughput_factor(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let temp_internal_air = 30.;
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: exhaust_runtime_throughput";

        let heat_source_wet_details = create_heat_pump_with_exhaust_input_from_json(2., 70);

        let mut heat_pump = create_heat_pump(
            heat_source_wet_details,
            energy_supply_conn_name_auxiliary,
            None,
            None,
            Some(101.),
            Some(temp_internal_air),
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        let (time_running, throughput_factor_zone) = heat_pump
            .running_time_throughput_factor(
                0.,
                "service_runtime_throughput",
                &HeatingServiceType::DomesticHotWaterRegular,
                1.,
                330.,
                330.,
                340.,
                design_flow_temp_op_cond_k(55.),
                1350.,
                true,
                100.,
                simulation_time_for_heat_pump.iter().current_iteration(),
                Some(TempSpreadCorrectionArg::Float(1.)),
                None,
            )
            .unwrap();

        assert_relative_eq!(time_running, 0.1305986895949177);
        assert_relative_eq!(throughput_factor_zone, 1.);
    }

    /// this test was added to guard against a deadlock issue with demo_hp_warm_air.json (use of temp_spread_correction_fn)
    #[rstest]
    fn test_demand_energy_on_heat_pump_service_space_warm_air(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: sink_air";
        let heat_pummp_sink_air = create_heat_pump_sink_air(
            energy_supply_conn_name_auxiliary,
            external_conditions,
            simulation_time_for_heat_pump,
        );
        let heat_pummp_sink_air = Arc::from(Mutex::from(heat_pummp_sink_air));

        let service_name = "service_space_warmair";
        let control = Arc::from(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(true)],
            0,
            1.,
        )));
        let volume_heated = 250.;
        let frac_convective = 0.9;

        // create_service_space_heating_warm_air expects a heat pump with sink air
        let mut heat_pump_service_space_warm_air = HeatPump::create_service_space_heating_warm_air(
            heat_pummp_sink_air.clone(),
            service_name,
            control,
            frac_convective,
            volume_heated,
        )
        .unwrap();

        let energy_demanded = 0.;

        let result = heat_pump_service_space_warm_air.demand_energy(
            energy_demanded,
            simulation_time_for_heat_pump.iter().current_iteration(),
        );

        assert!(result.is_ok());
    }

    /// this test was added to guard against a deadlock issue related to what we found with demo_hp_warm_air.json (use of temp_spread_correction_fn)
    #[rstest]
    fn test_running_time_throughput_factor_on_heat_pump_service_space(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: with_exhaust";
        let heat_pump_with_exhaust = create_heat_pump_with_exhaust(
            energy_supply_conn_name_auxiliary,
            None,
            None,
            external_conditions,
            simulation_time_for_heat_pump,
        );
        let heat_pump_with_exhaust = Arc::from(Mutex::from(heat_pump_with_exhaust));

        let service_name = "service_space";
        let temp_limit_upper = 50.0;
        let temp_diff_emit_design = 50.0;
        let control = Arc::from(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(true)],
            0,
            1.,
        )));
        let volume_heated = 250.0;

        let heat_pump_service_space = HeatPump::create_service_space_heating(
            heat_pump_with_exhaust.clone(),
            service_name,
            HeatPumpEmitterType::RadiatorsUfh,
            temp_limit_upper,
            temp_diff_emit_design,
            55.,
            control,
            volume_heated,
        );

        let space_heat_running_time_cumulative = 0.;
        let energy_demanded = 0.;
        let temp_flow = 0.;
        let temp_return = 0.;
        let result = heat_pump_service_space.running_time_throughput_factor(
            space_heat_running_time_cumulative,
            energy_demanded,
            temp_flow,
            temp_return,
            None,
            simulation_time_for_heat_pump.iter().current_iteration(),
        );

        assert!(result.is_ok());
    }

    // TODO: add more tests for other call sites of temp_spread_correction_fn:
    // HeatPumpServiceSpace: energy_output_max
    // HeatPumpServiceSpace: demand_energy
    // HeatPumpServiceSpaceWarmAir: running_time_throughput_factor (will this ever be reached though?)

    #[rstest]
    fn test_calc_throughput_factor(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        assert_eq!(heat_pump.calc_throughput_factor(10.), 1.0);
    }

    #[rstest]
    fn test_calc_energy_input(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let heat_pump = Arc::new(Mutex::new(create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        )));
        HeatPump::create_service_connection(heat_pump.clone(), "service_test").unwrap();

        let mut heat_pump = heat_pump.lock();

        let service_results = vec![
            ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
                service_name: "service_test".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 273.15,
                cop_op_cond: Some(2.9706605881515196),
                thermal_capacity_op_cond: Some(8.417674488123662),
                time_running_full_load: 0.11879765621857688,
                time_running_part_load: None,
                time_available_for_capacity_calc: 0.5,
                time_available_for_load_ratio_calc: 0.5,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 0.9999999999999999,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 0.9999999999999999,
                energy_heating_circ_pump: 0.001781964843278653,
                energy_source_circ_pump: 0.0011879765621857687,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            })),
            ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
                service_name: "service_test".into(),
                service_type: HeatingServiceType::Space,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 275.65,
                cop_op_cond: Some(3.091723311370327),
                thermal_capacity_op_cond: Some(6.9608039370972525),
                time_running_full_load: 0.07,
                time_running_part_load: None,
                time_available_for_capacity_calc: 0.38120234378142312,
                time_available_for_load_ratio_calc: 0.38120234378142312,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 1.,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 1.,
                energy_heating_circ_pump: 0.0022146353796754846,
                energy_source_circ_pump: 0.0014764235864503231,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            })),
            ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
                service_name: "service_test".into(),
                service_type: HeatingServiceType::Space,
                service_on: true,
                energy_output_required: 0.97,
                temp_output: Some(330.0),
                temp_source: 275.65,
                cop_op_cond: Some(3.091723311370327),
                thermal_capacity_op_cond: Some(6.9608039370972525),
                time_running_full_load: 0.05,
                time_running_part_load: None,
                #[allow(clippy::excessive_precision)]
                time_available_for_capacity_calc: 0.23355998513639082,
                #[allow(clippy::excessive_precision)]
                time_available_for_load_ratio_calc: 0.23355998513639082,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                time_constant_for_service: 1560.,
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 1.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 1.0,
                energy_heating_circ_pump: 0.002154923502450379,
                energy_source_circ_pump: 0.0014366156683002528,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            })),
        ];

        heat_pump.service_results = Arc::new(RwLock::new(service_results));

        let expected_results = vec![
            ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
                service_name: "service_test".into(),
                service_type: HeatingServiceType::DomesticHotWaterRegular,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 273.15,
                cop_op_cond: Some(2.9706605881515196),
                thermal_capacity_op_cond: Some(8.417674488123662),
                time_running_full_load: 0.11879765621857688,
                time_running_part_load: Some(0.11879765621857688),
                time_available_for_capacity_calc: 0.5,
                time_available_for_load_ratio_calc: 0.5,
                time_constant_for_service: 1560.,
                compressor_power_min_load: 2.83360358355901,
                load_ratio_continuous_min: 1.,
                load_ratio: 1.,
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: false,
                energy_input_hp: 0.33662546437937074,
                energy_delivered_hp: 0.9999999999999999,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: 0.33959540578483516,
                energy_delivered_total: 0.9999999999999999,
                energy_heating_circ_pump: 0.001781964843278653,
                energy_source_circ_pump: 0.0011879765621857687,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            })),
            ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
                service_name: "service_test".into(),
                service_type: HeatingServiceType::Space,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 275.65,
                cop_op_cond: Some(3.091723311370327),
                thermal_capacity_op_cond: Some(6.9608039370972525),
                time_running_full_load: 0.07,
                time_running_part_load: Some(0.17500000000000002),
                time_available_for_capacity_calc: 0.38120234378142315,
                time_available_for_load_ratio_calc: 0.38120234378142315,
                time_constant_for_service: 1560.,
                compressor_power_min_load: 0.9005726885711587,
                load_ratio_continuous_min: 0.4,
                load_ratio: 0.3147934475156495,
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: true,
                energy_input_hp: 0.13840863263212141,
                energy_delivered_hp: 1.,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: 0.1420996915982472,
                energy_delivered_total: 1.,
                energy_heating_circ_pump: 0.0022146353796754846,
                energy_source_circ_pump: 0.0014764235864503231,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            })),
            ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
                service_name: "service_test".into(),
                service_type: HeatingServiceType::Space,
                service_on: true,
                energy_output_required: 0.97,
                temp_output: Some(330.0),
                temp_source: 275.65,
                cop_op_cond: Some(3.091723311370327),
                thermal_capacity_op_cond: Some(6.9608039370972525),
                time_running_full_load: 0.05,
                time_running_part_load: Some(0.125),
                time_available_for_capacity_calc: 0.2335599851363908,
                time_available_for_load_ratio_calc: 0.2335599851363908,
                compressor_power_min_load: 0.9005726885711587,
                load_ratio_continuous_min: 0.4,
                load_ratio: 0.3147934475156495,
                time_constant_for_service: 1560.,
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: true,
                energy_input_hp: 0.13840863263212141,
                energy_delivered_hp: 1.0,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: 0.14200017180287205,
                energy_delivered_total: 1.0,
                energy_heating_circ_pump: 0.002154923502450379,
                energy_source_circ_pump: 0.0014366156683002528,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: None,
            })),
        ];

        heat_pump.calc_energy_input(0).unwrap();

        assert_eq!(*heat_pump.service_results.read().deref(), expected_results);
    }

    #[rstest]
    fn test_calculate_energy_input_error(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let mut heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );
        heat_pump.service_results = Arc::new(RwLock::new(vec![ServiceResult::Full(Box::new(
            HeatPumpEnergyCalculation {
                service_name: "service_water".into(),
                service_type: HeatingServiceType::DomesticHotWaterCombi,
                service_on: true,
                energy_output_required: 1.0,
                temp_output: Some(330.0),
                temp_source: 273.15,
                cop_op_cond: Some(2.9706605881515196),
                thermal_capacity_op_cond: Some(8.417674488123662),
                time_running_full_load: 0.11879765621857688,
                time_running_part_load: None,
                time_available_for_capacity_calc: 0.5,
                time_available_for_load_ratio_calc: 0.5,
                time_constant_for_service: 1560.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: false,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 0.9999999999999999,
                energy_input_backup: 0.0,
                energy_delivered_backup: 0.0,
                energy_input_total: Default::default(),
                energy_delivered_total: 0.9999999999999999,
                energy_heating_circ_pump: 0.001781964843278653,
                energy_source_circ_pump: 0.0011879765621857687,
                energy_output_required_boiler: 0.0,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: Some(0.),
            },
        ))]));

        assert!(heat_pump.calc_energy_input(0).is_err());
    }

    #[rstest]
    fn test_calc_auxiliary_energy(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let heat_pump = create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        );

        let (energy_standby, energy_crankcase_heater_mode, energy_off_mode) =
            heat_pump.calc_auxiliary_energy(1., 0.5, 0).unwrap();

        assert_eq!(energy_standby, 0.0);
        assert_eq!(energy_crankcase_heater_mode, 0.0);
        assert_eq!(energy_off_mode, 0.015);
    }

    #[rstest]
    fn test_extract_energy_from_source(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let energy_supply_conn_name_auxiliary = "HeatPump_auxiliary: hp1";

        let heat_pump_with_nw =
            Arc::new(Mutex::new(create_heat_pump_with_energy_supply_heat_source(
                energy_supply_conn_name_auxiliary,
                external_conditions,
                simulation_time_for_heat_pump,
            )));

        HeatPump::create_service_connection(heat_pump_with_nw.clone(), "service_extract_energy")
            .unwrap();

        let mut heat_pump_with_nw = heat_pump_with_nw.lock();

        let simtime = simulation_time_for_heat_pump.iter().current_iteration();

        heat_pump_with_nw
            .demand_energy(
                "service_extract_energy",
                &HeatingServiceType::DomesticHotWaterRegular,
                1.0,
                Some(330.0),
                330.0,
                340.0,
                design_flow_temp_op_cond_k(55.),
                1560.,
                true,
                simtime,
                Some(TempSpreadCorrectionArg::Float(1.0)),
                None,
                None,
                None,
                None,
                None,
            )
            .unwrap();

        let test_service_results = heat_pump_with_nw.service_results.read().clone();
        // Call the method under test
        heat_pump_with_nw.extract_energy_from_source(0).unwrap();
        let actual_service_results = heat_pump_with_nw.service_results.read().clone();
        assert_eq!(actual_service_results, test_service_results);

        // upstream uses mock to now check demand_energy is delegated to - was not considered that this is useful enough test to migrate

        // Test error when energy_delivered_hp < energy_input_hp
        if let Some(ServiceResult::Full(result)) =
            heat_pump_with_nw.service_results.write().get_mut(0)
        {
            result.energy_delivered_hp = 30.;
            result.energy_input_hp = 100.;
        }

        assert!(heat_pump_with_nw.extract_energy_from_source(0).is_err());
    }

    #[rstest]
    fn test_timestep_end(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        // Call demand_energy function to record the state of variables
        let heat_pump = Arc::new(Mutex::new(create_default_heat_pump(
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        )));

        HeatPump::create_service_connection(heat_pump.clone(), "servicetimestep_demand_energy")
            .unwrap();

        let mut heat_pump = heat_pump.lock();

        let simtime = simulation_time_for_heat_pump.iter().current_iteration();

        heat_pump
            .demand_energy(
                "servicetimestep_demand_energy",
                &HeatingServiceType::DomesticHotWaterRegular,
                5.0,
                Some(330.0),
                330.0,
                340.0,
                design_flow_temp_op_cond_k(55.),
                1560.,
                true,
                simtime,
                Some(TempSpreadCorrectionArg::Float(1.0)),
                None,
                None,
                None,
                None,
                None,
            )
            .unwrap();

        assert_relative_eq!(
            heat_pump.total_time_running_current_timestep_full_load,
            0.5939882810928845
        );

        let expected_service_results = [ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
            service_name: "servicetimestep_demand_energy".into(),
            service_type: HeatingServiceType::DomesticHotWaterRegular,
            service_on: true,
            energy_output_required: 5.0,
            temp_output: Some(330.0),
            temp_source: 273.15,
            cop_op_cond: Some(3.182442675928905), // 3.1824426759289044 in Python
            thermal_capacity_op_cond: Some(8.417674488123662),
            time_running_full_load: 0.5939882810928845,
            time_running_part_load: None,
            time_available_for_capacity_calc: 1.,
            time_available_for_load_ratio_calc: 1.,
            time_constant_for_service: 1560.,
            compressor_power_min_load: Default::default(),
            load_ratio_continuous_min: Default::default(),
            load_ratio: Default::default(),
            use_backup_heater_only: false,
            hp_operating_in_onoff_mode: false,
            energy_input_hp: Default::default(),
            energy_delivered_hp: 5.0,
            energy_input_backup: Default::default(),
            energy_delivered_backup: 0.0,
            energy_input_total: Default::default(),
            energy_delivered_total: 5.0,
            energy_heating_circ_pump: 0.008909824216393266,
            energy_source_circ_pump: 0.005939882810928845,
            energy_output_required_boiler: 0.0,
            energy_heating_warm_air_fan: 0.,
            energy_output_delivered_boiler: Some(0.0),
        }))];

        assert_eq!(
            *heat_pump.service_results.read().deref(),
            expected_service_results
        );

        // Call the method under test
        heat_pump.timestep_end(0).unwrap();

        assert_eq!(heat_pump.total_time_running_current_timestep_full_load, 0.);
        assert_eq!(*heat_pump.service_results.read().deref(), []);
    }

    // skipping Python's test_timestep_end_extract_energy due to mocking

    #[rstest]
    fn test_output_detailed_results_water(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
        energy_supply: EnergySupply,
    ) {
        let input = create_heat_pump_input_from_json(None, None);
        let heat_pump = Arc::new(Mutex::new(
            HeatPump::new(
                &input,
                Arc::from(RwLock::from(energy_supply)),
                "HeatPump_auxiliary: hp",
                simulation_time_for_heat_pump.step,
                Arc::new(external_conditions),
                2,
                None,
                None,
                true,
                None,
                None,
                create_temp_internal_air_fn(20.),
            )
            .unwrap(),
        ));

        HeatPump::create_service_connection(heat_pump.clone(), "servicetimestep_demand_energy")
            .unwrap();

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            heat_pump
                .lock()
                .demand_energy(
                    "servicetimestep_demand_energy",
                    &HeatingServiceType::DomesticHotWaterRegular,
                    5.,
                    Some(330.),
                    330.,
                    340.,
                    design_flow_temp_op_cond_k(55.),
                    1560.,
                    true,
                    t_it,
                    Some(TempSpreadCorrectionArg::Float(1.)),
                    None,
                    None,
                    None,
                    None,
                    None,
                )
                .unwrap();

            heat_pump.lock().timestep_end(t_idx).unwrap();
        }

        let mut expected_results_per_timestep: ResultsPerTimestep = indexmap! {
            "auxiliary".into() => indexmap! {
                ("energy_standby".into(), Some("kWh".into())) => vec![0.006090175783606732.into(), 0.006330405996936252.into()],
                ("energy_crankcase_heater_mode".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("energy_off_mode".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
            },
            "servicetimestep_demand_energy".into() => indexmap! {
                ("service_name".into(), None) => vec![ResultParamValue::String("servicetimestep_demand_energy".into()); 2],
                ("service_type".into(), None) => vec![ResultParamValue::String(HeatingServiceType::DomesticHotWaterRegular.to_string().into()); 2],
                ("service_on".into(), None) => vec![ResultParamValue::Boolean(true); 2],
                ("energy_output_required".into(), Some("kWh".into())) => vec![5.0.into(); 2],
                ("temp_output".into(), Some("K".into())) => vec![330.0.into(); 2],
                ("temp_source".into(), Some("K".into())) => vec![273.15.into(), 275.65.into()],
                ("thermal_capacity_op_cond".into(), Some("kW".into())) => vec![8.417674488123662.into(), 8.650924134797519.into()],
                ("cop_op_cond".into(), None) => vec![
                    3.182442675928905.into(), // 3.1824426759289044 in Python
                    3.1971341044162673.into(), // 3.197134104416267 in Python
                ],
                ("time_running_full_load".into(), Some("hours".into())) => vec![0.5939882810928845.into(), 0.5779729335375832.into()],
                ("time_running_part_load".into(), Some("hours".into())) => vec![0.5939882810928845.into(), 0.5779729335375832.into()],
                ("load_ratio".into(), None) => vec![1.0.into(); 2],
                ("hp_operating_in_onoff_mode".into(), None) => vec![ResultParamValue::Boolean(false); 2],
                ("energy_delivered_HP".into(), Some("kWh".into())) => vec![
                    5.0.into(),
                    5.000000000000001.into() // 5. in Python
                ],
                ("energy_delivered_backup".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("energy_delivered_total".into(), Some("kWh".into())) => vec![
                    5.0.into(),
                    5.000000000000001.into() // 5. in Python
                ],
                ("energy_input_HP".into(), Some("kWh".into())) => vec![
                    1.5711202083288363.into(), // 1.5711202083288365 in Python
                    1.5639006174603054.into()], // 1.5639006174603056 in Python
                ("energy_input_backup".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("energy_heating_circ_pump".into(), Some("kWh".into())) => vec![0.008909824216393266.into(), 0.008669594003063746.into()],
                ("energy_source_circ_pump".into(), Some("kWh".into())) => vec![0.005939882810928845.into(), 0.005779729335375832.into()],
                ("energy_heating_warm_air_fan".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("energy_input_total".into(), Some("kWh".into())) => vec![
                    1.5859699153561584.into(), // 1.5859699153561586 in Python
                    1.578349940798745.into() // 1.5783499407987451 in Python
                ],
                ("energy_output_delivered_boiler".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("energy_delivered_H5".into(), Some("kWh".into())) => vec![100.0.into()],
            },
        };

        let mut expected_results_annual: ResultsAnnual = indexmap! {
            "Overall".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 10.0.into(),
                ("time_running_full_load".into(), Some("hours".into())) => 1.1719612146304677.into(),
                ("time_running_part_load".into(), Some("hours".into())) => 1.1719612146304677.into(),
                ("energy_delivered_HP".into(), Some("kWh".into())) => 10.0.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.0.into(),
                ("energy_input_HP".into(), Some("kWh".into())) => 3.135020825789142.into(),
                ("energy_input_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_heating_circ_pump".into(), Some("kWh".into())) => 0.017579418219457014.into(),
                ("energy_source_circ_pump".into(), Some("kWh".into())) => 0.011719612146304677.into(),
                ("energy_heating_warm_air_fan".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_input_total".into(), Some("kWh".into())) => 3.1643198561549033.into(), // 3.1643198561549037 in Python
                ("energy_output_delivered_boiler".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_H5".into(), Some("kWh".into())) => 100.0.into(),
                ("CoP (H1)".into(), None) => 3.1771838471558898.into(),
                ("CoP (H2)".into(), None) => 3.1653973753129576.into(),
                ("CoP (H3)".into(), None) => 3.1653973753129576.into(),
                ("CoP (H4)".into(), None) => 3.1478807272334057.into(),
                ("CoP (H5)".into(), Some("Note: For water heating services, only valid when HP is only heat source".into())) => 31.478807272334056.into(),
            },
            "auxiliary".into() => indexmap! {
                ("energy_standby".into(), Some("kWh".into())) => 0.012420581780542984.into(),
                ("energy_crankcase_heater_mode".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_off_mode".into(), Some("kWh".into())) => 0.0.into(),
            },
            "servicetimestep_demand_energy".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 10.0.into(),
                ("time_running_full_load".into(), Some("hours".into())) => 1.1719612146304677.into(),
                ("time_running_part_load".into(), Some("hours".into())) => 1.1719612146304677.into(),
                ("energy_delivered_HP".into(), Some("kWh".into())) => 10.0.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.0.into(),
                ("energy_input_HP".into(), Some("kWh".into())) => 3.135020825789142.into(),
                ("energy_input_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_heating_circ_pump".into(), Some("kWh".into())) => 0.017579418219457014.into(),
                ("energy_source_circ_pump".into(), Some("kWh".into())) => 0.011719612146304677.into(),
                ("energy_heating_warm_air_fan".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_input_total".into(), Some("kWh".into())) => 3.1643198561549033.into(), // 3.1643198561549037 in Python
                ("energy_output_delivered_boiler".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_H5".into(), Some("kWh".into())) => 100.0.into(),
                ("CoP (H1)".into(), None) => 3.1897714738411085.into(),
                ("CoP (H2)".into(), None) => 3.1778915983807448.into(),
                ("CoP (H3)".into(), None) => 3.1778915983807448.into(),
                ("CoP (H4)".into(), None) => 3.160236782178972.into(),
                ("CoP (H5)".into(), Some("Note: For water heating services, only valid when HP is only heat source".into())) => 31.602367821789723.into(),
            },
        };

        let (results_per_timestep, results_annual) = heat_pump.lock().output_detailed_results(
            &indexmap! { "hwsname".into() => vec![100.0.into()] },
            &indexmap! { "servicetimestep_demand_energy".into() => "hwsname".into() },
        );

        assert_eq!(results_per_timestep, expected_results_per_timestep);
        assert_eq!(results_annual, expected_results_annual);

        // Test case where hot water source is not in hot water energy source data
        expected_results_per_timestep["servicetimestep_demand_energy"].insert(
            ("energy_delivered_H5".into(), Some("kWh".into())),
            vec![ResultParamValue::Empty; 2],
        );
        expected_results_annual["servicetimestep_demand_energy"].insert(
            ("energy_delivered_H5".into(), Some("kWh".into())),
            ResultParamValue::Empty,
        );
        expected_results_annual["Overall"].insert(
            ("energy_delivered_H5".into(), Some("kWh".into())),
            ResultParamValue::Empty,
        );
        expected_results_annual["servicetimestep_demand_energy"].shift_remove(&(
            "CoP (H5)".into(),
            Some("Note: For water heating services, only valid when HP is only heat source".into()),
        ));
        expected_results_annual["Overall"].shift_remove(&(
            "CoP (H5)".into(),
            Some("Note: For water heating services, only valid when HP is only heat source".into()),
        ));
        expected_results_annual["servicetimestep_demand_energy"].insert(
            (
                "CoP (H5)".into(),
                Some("Note: Cannot calculate CoP (H5) when HP is heating a pre-heat tank".into()),
            ),
            ResultParamValue::Empty,
        );
        expected_results_annual["Overall"].insert(
            (
                "CoP (H5)".into(),
                Some("Note: Cannot calculate CoP (H5) when HP is heating a pre-heat tank".into()),
            ),
            ResultParamValue::Empty,
        );

        let (results_per_timestep, results_annual) = heat_pump.lock().output_detailed_results(
            &indexmap! { "hwsname".into() => vec![100.0.into()] },
            &indexmap! { "servicetimestep_demand_energy".into() => "hwsname_other".into()},
        );

        assert_eq!(results_per_timestep, expected_results_per_timestep);
        assert_eq!(results_annual, expected_results_annual);
    }

    #[rstest]
    fn test_output_detailed_results_space(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
        energy_supply: EnergySupply,
    ) {
        let input = create_heat_pump_input_from_json(None, None);
        let heat_pump = Arc::new(Mutex::new(
            HeatPump::new(
                &input,
                Arc::from(RwLock::from(energy_supply)),
                "HeatPump_auxiliary: hp",
                simulation_time_for_heat_pump.step,
                Arc::new(external_conditions),
                2,
                None,
                None,
                true,
                None,
                None,
                create_temp_internal_air_fn(20.),
            )
            .unwrap(),
        ));

        HeatPump::create_service_connection(heat_pump.clone(), "servicetimestep_demand_energy")
            .unwrap();

        for (t_idx, t_it) in simulation_time_for_heat_pump.iter().enumerate() {
            heat_pump
                .lock()
                .demand_energy(
                    "servicetimestep_demand_energy",
                    &HeatingServiceType::Space,
                    5.,
                    Some(330.),
                    330.,
                    340.,
                    design_flow_temp_op_cond_k(55.),
                    1560.,
                    true,
                    t_it,
                    Some(TempSpreadCorrectionArg::Float(1.)),
                    None,
                    None,
                    None,
                    None,
                    None,
                )
                .unwrap();

            heat_pump.lock().timestep_end(t_idx).unwrap();
        }

        let expected_results_per_timestep: ResultsPerTimestep = indexmap! {
            "auxiliary".into() => indexmap! {
                ("energy_standby".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("energy_crankcase_heater_mode".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("energy_off_mode".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
            },
            "servicetimestep_demand_energy".into() => indexmap! {
                ("service_name".into(), None) => vec![ResultParamValue::String("servicetimestep_demand_energy".into()); 2],
                ("service_type".into(), None) => vec![ResultParamValue::String(HeatingServiceType::Space.to_string().into()); 2],
                ("service_on".into(), None) => vec![ResultParamValue::Boolean(true); 2],
                ("energy_output_required".into(), Some("kWh".into())) => vec![5.0.into(); 2],
                ("temp_output".into(), Some("K".into())) => vec![330.0.into(); 2],
                ("temp_source".into(), Some("K".into())) => vec![273.15.into(), 275.65.into()],
                ("thermal_capacity_op_cond".into(), Some("kW".into())) => vec![8.417674488123662.into(), 8.650924134797519.into()],
                ("cop_op_cond".into(), None) => vec![
                    3.182442675928905.into(), // 3.1824426759289044 in Python
                    3.1971341044162673.into(), // 3.197134104416267 in Python
                ],
                ("time_running_full_load".into(), Some("hours".into())) => vec![0.5939882810928845.into(), 0.5779729335375832.into()],
                ("time_running_part_load".into(), Some("hours".into())) => vec![1.0.into(), 1.0.into()],
                ("load_ratio".into(), None) => vec![0.5939882810928845.into(), 0.5779729335375832.into()],
                ("hp_operating_in_onoff_mode".into(), None) => vec![ResultParamValue::Boolean(false); 2],
                ("energy_delivered_HP".into(), Some("kWh".into())) => vec![
                    5.0.into(),
                    5.000000000000001.into()
                ],
                ("energy_delivered_backup".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("energy_delivered_total".into(), Some("kWh".into())) => vec![
                    5.0.into(),
                    5.000000000000001.into()
                ],
                ("energy_input_HP".into(), Some("kWh".into())) => vec![
                    1.5711202083288363.into(), // 1.5711202083288365 in Python
                    1.5639006174603054.into()], // 1.5639006174603056 in Python
                ("energy_input_backup".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("energy_heating_circ_pump".into(), Some("kWh".into())) => vec![0.008909824216393266.into(), 0.008669594003063746.into()],
                ("energy_source_circ_pump".into(), Some("kWh".into())) => vec![0.005939882810928845.into(), 0.005779729335375832.into()],
                ("energy_heating_warm_air_fan".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("energy_input_total".into(), Some("kWh".into())) => vec![
                    1.5859699153561584.into(), // 1.5859699153561586 in Python
                    1.578349940798745.into() // 1.5783499407987451 in Python
                ],
                ("energy_output_delivered_boiler".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("energy_delivered_H5".into(), Some("kWh".into())) => vec![5.0.into(), 5.000000000000001.into()],
            },
        };

        let expected_results_annual: ResultsAnnual = indexmap! {
            "Overall".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 10.0.into(),
                ("time_running_full_load".into(), Some("hours".into())) => 1.1719612146304677.into(),
                ("time_running_part_load".into(), Some("hours".into())) => 2.0.into(),
                ("energy_delivered_HP".into(), Some("kWh".into())) => 10.0.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.0.into(),
                ("energy_input_HP".into(), Some("kWh".into())) => 3.135020825789142.into(),
                ("energy_input_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_heating_circ_pump".into(), Some("kWh".into())) => 0.017579418219457014.into(),
                ("energy_source_circ_pump".into(), Some("kWh".into())) => 0.011719612146304677.into(),
                ("energy_heating_warm_air_fan".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_input_total".into(), Some("kWh".into())) => 3.1643198561549033.into(), // 3.1643198561549037 in Python
                ("energy_output_delivered_boiler".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_H5".into(), Some("kWh".into())) => 10.0.into(),
                ("CoP (H1)".into(), None) => 3.1897714738411085.into(),
                ("CoP (H2)".into(), None) => 3.1778915983807448.into(),
                ("CoP (H3)".into(), None) => 3.1778915983807448.into(),
                ("CoP (H4)".into(), None) => 3.160236782178972.into(),
                ("CoP (H5)".into(), Some("Note: For water heating services, only valid when HP is only heat source".into())) => 3.160236782178972.into(),
            },
            "auxiliary".into() => indexmap! {
                ("energy_standby".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_crankcase_heater_mode".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_off_mode".into(), Some("kWh".into())) => 0.0.into(),
            },
            "servicetimestep_demand_energy".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 10.0.into(),
                ("time_running_full_load".into(), Some("hours".into())) => 1.1719612146304677.into(),
                ("time_running_part_load".into(), Some("hours".into())) => 2.0.into(),
                ("energy_delivered_HP".into(), Some("kWh".into())) => 10.0.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.0.into(),
                ("energy_input_HP".into(), Some("kWh".into())) => 3.135020825789142.into(),
                ("energy_input_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_heating_circ_pump".into(), Some("kWh".into())) => 0.017579418219457014.into(),
                ("energy_source_circ_pump".into(), Some("kWh".into())) => 0.011719612146304677.into(),
                ("energy_heating_warm_air_fan".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_input_total".into(), Some("kWh".into())) => 3.1643198561549033.into(), // 3.1643198561549037 in Python
                ("energy_output_delivered_boiler".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_H5".into(), Some("kWh".into())) => 10.0.into(),
                ("CoP (H1)".into(), None) => 3.1897714738411085.into(),
                ("CoP (H2)".into(), None) => 3.1778915983807448.into(),
                ("CoP (H3)".into(), None) => 3.1778915983807448.into(),
                ("CoP (H4)".into(), None) => 3.160236782178972.into(),
                ("CoP (H5)".into(), Some("Note: For water heating services, only valid when HP is only heat source".into())) => 3.160236782178972.into(),
            },
        };

        let (results_per_timestep, results_annual) = heat_pump.lock().output_detailed_results(
            &indexmap! { "hwsname".into() => vec![100.0.into()] },
            &indexmap! {},
        );

        assert_eq!(results_per_timestep, expected_results_per_timestep);
        assert_eq!(results_annual, expected_results_annual);
    }

    #[rstest]
    fn test_calc_service_cop_zero_values(
        external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
        energy_supply: EnergySupply,
    ) {
        let input = create_heat_pump_input_from_json(None, None);
        let heat_pump = Arc::new(Mutex::new(
            HeatPump::new(
                &input,
                Arc::from(RwLock::from(energy_supply)),
                "HeatPump_auxiliary: hp",
                simulation_time_for_heat_pump.step,
                Arc::new(external_conditions),
                2,
                None,
                None,
                true,
                None,
                None,
                create_temp_internal_air_fn(20.),
            )
            .unwrap(),
        ));

        let mut results_totals: ResultAnnual = indexmap! {
            ("energy_delivered_HP".into(), Some("kWh".into())) => 0.0.into(),
            ("energy_input_HP".into(), Some("kWh".into())) => 0.0.into(),
            ("energy_source_circ_pump".into(), Some("kWh".into())) => 0.0.into(),
            ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
            ("energy_input_backup".into(), Some("kWh".into())) => 0.0.into(),
            ("energy_delivered_H5".into(), Some("kWh".into())) => 0.0.into(),
            ("energy_heating_circ_pump".into(), Some("kWh".into())) => 0.0.into(),
            ("energy_heating_warm_air_fan".into(), Some("kWh".into())) => 0.0.into(),
        };

        heat_pump.lock().calc_service_cop(&mut results_totals, None);

        assert_eq!(results_totals[&("CoP (H1)".into(), None)], 0.);
        assert_eq!(results_totals[&("CoP (H2)".into(), None)], 0.);
        assert_eq!(results_totals[&("CoP (H3)".into(), None)], 0.);
        assert_eq!(results_totals[&("CoP (H4)".into(), None)], 0.);
        assert_eq!(results_totals[&("CoP (H5)".into(), None)], 0.);
    }

    #[rstest]
    fn test_backup_only_operation(
        mut external_conditions: ExternalConditions,
        simulation_time_for_heat_pump: SimulationTime,
    ) {
        let temp_limit_upper = 65.;
        let control = Arc::new(create_setpoint_time_control(vec![Some(20.), Some(20.)]));
        external_conditions.air_temps = vec![-10., 2.5];
        let heat_pump_input = create_heat_pump_input_from_json(Some("Substitute"), None);
        let heat_pump = Arc::new(Mutex::new(create_heat_pump(
            heat_pump_input,
            "hp_aux",
            None,
            None,
            Some(101.),
            None,
            external_conditions,
            simulation_time_for_heat_pump,
            None,
        )));
        let hp_service_space1 = HeatPump::create_service_space_heating(
            heat_pump.clone(),
            "hp_space_heating_1",
            HeatPumpEmitterType::RadiatorsUfh,
            temp_limit_upper,
            5.,
            55.,
            control.clone(),
            120.,
        );
        let hp_service_space2 = HeatPump::create_service_space_heating(
            heat_pump.clone(),
            "hp_space_heating_2",
            HeatPumpEmitterType::RadiatorsUfh,
            temp_limit_upper,
            5.,
            55.,
            control,
            120.,
        );
        hp_service_space1
            .demand_energy(
                2.,
                35.,
                30.,
                Some(0.),
                None,
                None,
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();
        hp_service_space2
            .demand_energy(
                2.,
                35.,
                30.,
                Some(0.25),
                None,
                None,
                simulation_time_for_heat_pump.iter().current_iteration(),
            )
            .unwrap();

        let expected_results = vec![
            ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
                service_name: "hp_space_heating_1".into(),
                service_type: HeatingServiceType::Space,
                service_on: true,
                energy_output_required: 2.0,
                temp_output: Some(308.15),
                temp_source: 263.15,
                cop_op_cond: Some(3.317191947724431), // 3.3171919477244307 in Python
                thermal_capacity_op_cond: Some(9.243872347406562),
                time_running_full_load: 0.6666666666666666,
                time_running_part_load: None,
                time_available_for_capacity_calc: 1.,
                time_available_for_load_ratio_calc: 1.,
                time_constant_for_service: 1370.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: true,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 0.,
                energy_input_backup: 2.,
                energy_delivered_backup: 2.,
                energy_input_total: Default::default(),
                energy_delivered_total: 2.,
                energy_heating_circ_pump: 0.009999999999999998,
                energy_source_circ_pump: 0.,
                energy_output_required_boiler: 0.,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: Some(0.),
            })),
            ServiceResult::Full(Box::new(HeatPumpEnergyCalculation {
                service_name: "hp_space_heating_2".into(),
                service_type: HeatingServiceType::Space,
                service_on: true,
                energy_output_required: 2.0,
                temp_output: Some(308.15),
                temp_source: 263.15,
                cop_op_cond: Some(3.317191947724431), // 3.3171919477244307 in Python
                thermal_capacity_op_cond: Some(9.243872347406562),
                time_running_full_load: 0.25,
                time_running_part_load: None,
                time_available_for_capacity_calc: 0.25,
                time_available_for_load_ratio_calc: 0.33333333333333337,
                time_constant_for_service: 1370.,
                compressor_power_min_load: Default::default(),
                load_ratio_continuous_min: Default::default(),
                load_ratio: Default::default(),
                use_backup_heater_only: true,
                hp_operating_in_onoff_mode: Default::default(),
                energy_input_hp: Default::default(),
                energy_delivered_hp: 0.,
                energy_input_backup: 0.75,
                energy_delivered_backup: 0.75,
                energy_input_total: Default::default(),
                energy_delivered_total: 0.75,
                energy_heating_circ_pump: 0.00375,
                energy_source_circ_pump: 0.,
                energy_output_required_boiler: 0.,
                energy_heating_warm_air_fan: 0.,
                energy_output_delivered_boiler: Some(0.),
            })),
        ];

        assert_eq!(
            *heat_pump.lock().service_results.read().deref(),
            expected_results
        );
    }

    // In Python below tests are in a separate class called TestHeatPump_HWOnly
    fn create_heat_pump_hw_only(
        vol_daily_average: Option<f64>,
        test_data: Option<HeatPumpHotWaterTestData>,
        simulation_time_for_heat_pump: SimulationTime,
    ) -> HeatPumpHotWaterOnly {
        let mut energy_supply = energy_supply(simulation_time_for_heat_pump);
        let test_data: HeatPumpHotWaterTestData = test_data.unwrap_or(HeatPumpHotWaterTestData {
            l: Some(HeatPumpHotWaterOnlyTestDatum {
                cop_dhw: 2.5,
                hw_tapping_prof_daily_total: 11.655,
                energy_input_measured: 4.6,
                power_standby: 0.03,
                hw_vessel_loss_daily: 1.6,
            }),
            m: HeatPumpHotWaterOnlyTestDatum {
                cop_dhw: 2.7,
                hw_tapping_prof_daily_total: 5.845,
                energy_input_measured: 2.15,
                power_standby: 0.02,
                hw_vessel_loss_daily: 1.18,
            },
        });

        let power_max = 3.0;
        let vol_daily_average = vol_daily_average.unwrap_or(150.);
        let tank_volume = 200.0;
        let daily_losses = 1.5;
        let heat_exchanger_surface_area = 1.2;
        let in_use_factor_mismatch = 0.6;
        let tank_volume_declared = 180.0;
        let heat_exchanger_surface_area_declared = 1.0;
        let daily_losses_declared = 1.2;

        energy_supply.register_end_user_name("end_user_name".into());
        let energy_supply_connection = EnergySupplyConnection::new(
            Arc::from(RwLock::from(energy_supply)),
            "end_user_name".into(),
        );
        let control_min = Arc::new(create_setpoint_time_control(vec![
            Some(52.),
            Some(52.),
            None,
            Some(52.),
        ]));
        let control_max = Arc::new(create_setpoint_time_control(vec![
            Some(60.),
            Some(60.),
            Some(60.),
            Some(60.),
        ]));

        HeatPumpHotWaterOnly::new(
            power_max,
            energy_supply_connection,
            &test_data,
            vol_daily_average,
            tank_volume,
            daily_losses,
            heat_exchanger_surface_area,
            in_use_factor_mismatch,
            tank_volume_declared,
            heat_exchanger_surface_area_declared,
            daily_losses_declared,
            1.,
            control_min,
            control_max,
        )
    }

    #[rstest]
    fn test_init_efficiencies(simulation_time_for_heat_pump: SimulationTime) {
        let test_data: HeatPumpHotWaterTestData = HeatPumpHotWaterTestData {
            l: None,
            m: HeatPumpHotWaterOnlyTestDatum {
                cop_dhw: 2.7,
                hw_tapping_prof_daily_total: 5.845,
                energy_input_measured: 2.15,
                power_standby: 0.02,
                hw_vessel_loss_daily: 1.18,
            },
        };
        let heat_pump =
            create_heat_pump_hw_only(None, Some(test_data), simulation_time_for_heat_pump);

        assert_eq!(heat_pump.initial_efficiency, 3.0478653375963884);

        let heat_pump = create_heat_pump_hw_only(Some(90.), None, simulation_time_for_heat_pump);

        assert_eq!(heat_pump.initial_efficiency, 3.0478653375963884);

        let heat_pump = create_heat_pump_hw_only(Some(200.), None, simulation_time_for_heat_pump);

        assert_eq!(heat_pump.initial_efficiency, 2.7473226825842696);
    }

    // skipping Python's test_init_efficiencies_invalid as test_data can't be empty in Rust

    #[rstest]
    fn test_calc_efficiency(simulation_time_for_heat_pump: SimulationTime) {
        let heat_pump = create_heat_pump_hw_only(None, None, simulation_time_for_heat_pump);

        assert_relative_eq!(
            heat_pump.calc_efficiency(),
            1.738556406,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_calc_efficiency_criteria(simulation_time_for_heat_pump: SimulationTime) {
        let mut heat_pump = create_heat_pump_hw_only(None, None, simulation_time_for_heat_pump);
        heat_pump.daily_losses = 1.1;

        assert_relative_eq!(heat_pump.calc_efficiency(), 2.8975940100903292);
    }

    #[rstest]
    fn test_setpnt(simulation_time_for_heat_pump: SimulationTime) {
        let simtime = simulation_time_for_heat_pump.iter().current_iteration();
        let heat_pump = create_heat_pump_hw_only(None, None, simulation_time_for_heat_pump);

        let (minsetpnt, maxsetpnt) = heat_pump.setpnt(simtime).unwrap();

        assert_eq!(minsetpnt, Some(52.));
        assert_eq!(maxsetpnt, Some(60.));
    }

    #[rstest]
    fn test_setpnt_errors(simulation_time_for_heat_pump: SimulationTime) {
        let simtime = simulation_time_for_heat_pump.iter().current_iteration();
        let heat_pump = create_heat_pump_hw_only(None, None, simulation_time_for_heat_pump);

        let mut hp1 = heat_pump.clone();
        hp1.control_min = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(true)],
            0,
            1.,
        )));

        assert!(hp1.setpnt(simtime).is_err());

        let mut hp2 = heat_pump.clone();
        hp2.control_max = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(true)],
            0,
            1.,
        )));

        assert!(hp2.setpnt(simtime).is_err());
    }

    #[rstest]
    fn test_demand_energy_for_hw_only(simulation_time_for_heat_pump: SimulationTime) {
        let simtime = simulation_time_for_heat_pump.iter().current_iteration();
        let heat_pump = create_heat_pump_hw_only(None, None, simulation_time_for_heat_pump);

        assert_relative_eq!(heat_pump.demand_energy(10., 50., 40., simtime), 3.);
    }

    #[rstest]
    fn test_demand_energy_off_for_hw_only(simulation_time_for_heat_pump: SimulationTime) {
        let simtime = simulation_time_for_heat_pump.iter().current_iteration();
        let mut heat_pump = create_heat_pump_hw_only(None, None, simulation_time_for_heat_pump);
        heat_pump.control_min = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(false)],
            0,
            1.,
        )));

        assert_relative_eq!(heat_pump.demand_energy(10., 50., 40., simtime), 0.);
    }

    #[rstest]
    fn test_energy_output_max_for_hw_only(simulation_time_for_heat_pump: SimulationTime) {
        let simtime = simulation_time_for_heat_pump.iter().current_iteration();
        let heat_pump = create_heat_pump_hw_only(None, None, simulation_time_for_heat_pump);

        assert_relative_eq!(heat_pump.energy_output_max(50., simtime), 3.);
    }

    #[rstest]
    fn test_energy_output_max_off(simulation_time_for_heat_pump: SimulationTime) {
        let simtime = simulation_time_for_heat_pump.iter().current_iteration();
        let mut heat_pump = create_heat_pump_hw_only(None, None, simulation_time_for_heat_pump);

        heat_pump.control_min = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(false)],
            0,
            1.,
        )));

        assert_relative_eq!(heat_pump.energy_output_max(50., simtime), 0.);
    }
}
