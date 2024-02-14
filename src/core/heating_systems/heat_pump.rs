use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::units::{celsius_to_kelvin, kelvin_to_celsius, HOURS_PER_DAY};
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::external_conditions::ExternalConditions;
use crate::input::{
    HeatPumpBackupControlType, HeatPumpHotWaterOnlyTestDatum, HeatPumpHotWaterTestData,
    HeatPumpSinkType, HeatPumpSourceType, HeatPumpTestDatum, HeatSourceWetDetails, TestLetter,
};
use crate::simulation_time::SimulationTimeIteration;
use arrayvec::ArrayString;
use derivative::Derivative;
use interp::interp;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use polyfit_rs::polyfit_rs::polyfit;
use serde::Deserialize;
use serde_enum_str::Serialize_enum_str;
use std::collections::HashMap;
use std::hash::Hash;
use std::iter::Sum;
use std::ops::{Add, Div};
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

/// This module provides objects to represent heat pumps and heat pump test data.
/// The calculations are based on the DAHPSE method developed for generating PCDB
/// entries for SAP 2012 and SAP 10. DAHPSE was based on a draft of
/// BS EN 15316-4-2:2017 and is described in the SAP calculation method CALCM-01.

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

#[derive(Copy, Clone, Serialize_enum_str)]
enum ServiceType {
    Water,
    Space,
}

/// Calculate Carnot CoP based on source and outlet temperatures (in Kelvin)
fn carnot_cop(temp_source: f64, temp_outlet: f64, temp_diff_limit_low: Option<f64>) -> f64 {
    let mut temp_diff = temp_outlet - temp_source;
    if let Some(low_limit) = temp_diff_limit_low {
        temp_diff = *[temp_diff, low_limit]
            .iter()
            .max_by(|a, b| a.total_cmp(b))
            .unwrap();
    }
    temp_outlet / temp_diff
}

/// Interpolate between test data records for different air flow rates
///    
/// Arguments:
/// * `throughput_exhaust_air` - throughput (litres / second) of exhaust air
/// * `test_data`
///        - list of records of heat pump test data, each with the following elements:
///                - air_flow_rate
///                - test_letter
///                - capacity
///                - cop
///                - degradation_coeff
///                - design_flow_temp (in Celsius)
///                - temp_outlet (in Celsius)
///                - temp_source (in Celsius)
///                - temp_test (in Celsius)
fn interpolate_exhaust_air_heat_pump_test_data(
    throughput_exhaust_air: f64,
    test_data: &Vec<HeatPumpTestDatum>,
) -> Result<(f64, Vec<HeatPumpTestDatum>), String> {
    // split test records into different lists by air flow rate
    let mut test_data_by_air_flow_rate: HashMap<OrderedFloat<f64>, Vec<&HeatPumpTestDatum>> =
        Default::default();
    for test_data_record in test_data {
        if test_data_record.air_flow_rate.is_none() {
            continue;
        }
        test_data_by_air_flow_rate
            .entry(OrderedFloat(test_data_record.air_flow_rate.unwrap()))
            .or_default()
            .push(test_data_record);
    }

    // check that all lists have same combination of design flow temp and test letter
    let mut fixed_temps_and_test_letters: Option<Vec<TestDatumTempsAndTestLetters>> = None;
    for (air_flow_rate, test_data_record_list) in &test_data_by_air_flow_rate {
        // find and save all the combinations of design flow temp and test letter for this air flow rate
        let mut fixed_temps_and_test_letters_this: Vec<TestDatumTempsAndTestLetters> =
            Default::default();
        for test_data_record in test_data_record_list {
            fixed_temps_and_test_letters_this.push(TestDatumTempsAndTestLetters {
                air_flow_rate: air_flow_rate.0,
                design_flow_temp: test_data_record.design_flow_temp,
                test_letter: &test_data_record.test_letter,
                temp_outlet: test_data_record.temp_outlet,
                temp_source: test_data_record.temp_source,
                temp_test: test_data_record.temp_test,
            })
        }

        match fixed_temps_and_test_letters {
            None => {
                fixed_temps_and_test_letters = Some(fixed_temps_and_test_letters_this);
            }
            Some(ref fixed_temps_and_test_letters) => {
                if (fixed_temps_and_test_letters != &fixed_temps_and_test_letters_this) {
                    return Err("In heat pump test data, fixed temps and test_letters were not consistent for one air_flow_temp value".to_string());
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
        .ok_or("Non-empty test data was expected".to_string())?
        .iter()
        .map(|datum| {
            // create lists of test data values ordered by air flow rate
            let mut capacity_list: Vec<f64> = Default::default();
            let mut cop_list: Vec<f64> = Default::default();
            let mut degradation_coeff_list: Vec<f64> = Default::default();
            for air_flow_rate in air_flow_rates_ordered {
                for test_record in &test_data_by_air_flow_rate[&OrderedFloat(*air_flow_rate)] {
                    if test_record.design_flow_temp == datum.design_flow_temp
                        && test_record.test_letter.as_str() == datum.test_letter
                    {
                        capacity_list.push(test_record.capacity);
                        cop_list.push(test_record.cop);
                        degradation_coeff_list.push(test_record.degradation_coefficient);
                    }
                }
            }

            let capacity = interp(
                air_flow_rates_ordered,
                &capacity_list,
                throughput_exhaust_air,
            );
            let cop = interp(air_flow_rates_ordered, &cop_list, throughput_exhaust_air);
            let degradation_coeff = interp(
                air_flow_rates_ordered,
                &degradation_coeff_list,
                throughput_exhaust_air,
            );

            HeatPumpTestDatum {
                air_flow_rate: Some(datum.air_flow_rate),
                test_letter: test_letter(datum.test_letter),
                capacity,
                cop,
                degradation_coefficient: degradation_coeff,
                design_flow_temp: datum.design_flow_temp,
                temp_outlet: datum.temp_outlet,
                temp_source: datum.temp_source,
                temp_test: datum.temp_test,
            }
        })
        .collect::<Vec<_>>();

    let lowest_air_flow_rate_in_test_data = &air_flow_rates_ordered
        .iter()
        .max_by(|a, b| a.total_cmp(b).reverse())
        .ok_or("Non-empty test data was expected".to_string())?;

    Ok((
        lowest_air_flow_rate_in_test_data.to_owned().to_owned(),
        test_data_interp_by_air_flow_rate,
    ))
}

fn test_letter(letter: &str) -> TestLetter {
    let mut test_letter = TestLetter::new();
    test_letter.push_str(letter);
    test_letter
}

#[derive(Derivative)]
#[derivative(PartialEq, Debug)]
struct TestDatumTempsAndTestLetters<'a> {
    #[derivative(PartialEq = "ignore")]
    // we need to ignore air flow rate when performing comparisons
    pub air_flow_rate: f64,
    pub design_flow_temp: f64,
    pub test_letter: &'a str,
    pub temp_outlet: f64,
    pub temp_source: f64,
    pub temp_test: f64,
}

#[derive(Copy, Clone, Debug, PartialEq)]
struct CompleteHeatPumpTestDatum {
    pub air_flow_rate: Option<f64>,
    pub test_letter: TestLetter,
    pub capacity: f64,
    pub cop: f64,
    pub degradation_coefficient: f64,
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

enum DatumItem {
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
            degradation_coefficient,
            design_flow_temp,
            temp_outlet,
            temp_source,
            temp_test,
        } = self;
        CompleteHeatPumpTestDatum {
            air_flow_rate: *air_flow_rate,
            test_letter: *test_letter,
            capacity: *capacity,
            cop: *cop,
            degradation_coefficient: *degradation_coefficient,
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
#[derive(Clone)]
struct HeatPumpTestData {
    test_data: HashMap<OrderedFloat<f64>, Vec<CompleteHeatPumpTestDatum>>,
    dsgn_flow_temps: Vec<OrderedFloat<f64>>,
    average_deg_coeff: Vec<f64>,
    average_cap: Vec<f64>,
    temp_spread_test_conditions: Vec<f64>,
    regression_coeffs: HashMap<OrderedFloat<f64>, Vec<f64>>,
}

const TEST_LETTERS_NON_BIVALENT: [char; 4] = ['A', 'B', 'C', 'D'];
const TEST_LETTERS_ALL: [char; 5] = ['A', 'B', 'C', 'D', 'F'];

impl HeatPumpTestData {
    pub fn new(data: Vec<HeatPumpTestDatum>) -> Result<Self, String> {
        // keyed by design flow temp
        let mut test_data: HashMap<OrderedFloat<f64>, Vec<HeatPumpTestDatum>> = Default::default();
        let mut dsgn_flow_temps: Vec<OrderedFloat<f64>> = Default::default();

        // variable to count duplicate records for each design flow temp
        let mut dupl: HashMap<OrderedFloat<f64>, usize> = Default::default();

        for datum in data {
            let mut saved_datum = datum.clone();
            let dsgn_flow_temp = OrderedFloat(datum.design_flow_temp);

            // When a new design flow temp is encountered, add it to lists/ maps
            if !dupl.contains_key(&dsgn_flow_temp) {
                if !dsgn_flow_temps.contains(&dsgn_flow_temp) {
                    dsgn_flow_temps.push(dsgn_flow_temp);
                }
                test_data
                    .entry(dsgn_flow_temp)
                    .or_insert_with(|| Default::default());
            }

            let mut duplicate = false;
            for d in test_data.get(&dsgn_flow_temp).unwrap() {
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
                *dupl.entry(dsgn_flow_temp).or_default() += 1;
            }

            test_data
                .entry(dsgn_flow_temp)
                .or_default()
                .push(saved_datum);
        }

        // Check the number of test records is as expected
        // - 1 or 2 design flow temps
        // - 4 or 5 distinct records for each flow temp
        if dsgn_flow_temps.is_empty() {
            return Err("No test data provided for heat pump performance".into());
        } else if dsgn_flow_temps.len() > 2 {
            return Err(
                "Test data for a maximum of 2 design flow temperatures may be provided".into(),
            );
        }
        for (dsgn_flow_temp, data) in test_data.iter() {
            if dupl.get(dsgn_flow_temp).is_some() {
                if (data.len() - dupl.get(dsgn_flow_temp).unwrap()) != 4 {
                    return Err(
                        "Expected 4 distinct records for each design flow temperature".into(),
                    );
                }
            } else if data.len() != 5 {
                return Err("Expected 5 records for each design flow temperature".into());
            }
        }

        // Check if test letters ABCDE are present as expected
        let mut test_letter_vec: Vec<char> = Default::default();
        for temperature in &dsgn_flow_temps {
            for test_data in &test_data[temperature] {
                for test_letter in test_data.test_letter.chars() {
                    test_letter_vec.push(test_letter);
                }
                if test_letter_vec.len() == 5 {
                    for test_letter_check in TEST_LETTERS_ALL {
                        if !test_letter_vec.contains(&test_letter_check) {
                            return Err(format!(
                                "Expected test letter {test_letter_check} in {temperature} degree temp data"
                            ));
                        }
                    }
                    test_letter_vec = Default::default();
                }
            }
        }

        dsgn_flow_temps.sort_by(|a, b| a.total_cmp(b));
        for (_, data) in test_data.iter_mut() {
            data.sort_by(|a, b| a.temp_test.total_cmp(&b.temp_test));
        }

        let average_deg_coeff = ave_degradation_coeff(&dsgn_flow_temps, &test_data);
        let average_cap = ave_capacity(&dsgn_flow_temps, &test_data);
        let temp_spread_test_conditions = init_temp_spread_test_conditions(&dsgn_flow_temps)?;
        let regression_coeffs = init_regression_coeffs(&dsgn_flow_temps, &test_data)?;

        // Calculate derived variables for each data record which are not time-dependent
        let test_data: HashMap<OrderedFloat<f64>, Vec<CompleteHeatPumpTestDatum>> = test_data
            .iter()
            .map(|(dsgn_flow_temp, data)| {
                let (carnot_cops, exergetic_effs) =
                    data.iter().fold((vec![], vec![]), |mut acc, datum| {
                        // Get the source and outlet temperatures from the test record
                        let temp_source = celsius_to_kelvin(datum.temp_source);
                        let temp_outlet = celsius_to_kelvin(datum.temp_outlet);

                        // Calculate the Carnot CoP
                        let carnot_cop = carnot_cop(temp_source, temp_outlet, None);
                        acc.0.push(carnot_cop);
                        acc.1.push(datum.cop / carnot_cop);
                        acc
                    });

                let temp_source_cld = celsius_to_kelvin(data[0].temp_source);
                let temp_outlet_cld = celsius_to_kelvin(data[0].temp_outlet);
                let carnot_cop_cld = carnot_cops[0];

                let theoretical_load_ratios =
                    data.iter().enumerate().fold(vec![], |mut acc, (i, datum)| {
                        // Get the source and outlet temperatures from the test record
                        let temp_source = celsius_to_kelvin(datum.temp_source);
                        let temp_outlet = celsius_to_kelvin(datum.temp_outlet);

                        let theoretical_load_ratio = ((carnot_cops[i] / carnot_cop_cld)
                            * (temp_outlet_cld * temp_source / (temp_source_cld * temp_outlet))
                                .powf(N_EXER));
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
                (*dsgn_flow_temp, complete_data)
            })
            .collect();

        Ok(Self {
            test_data,
            dsgn_flow_temps,
            average_deg_coeff,
            average_cap,
            temp_spread_test_conditions,
            regression_coeffs,
        })
    }

    /// Return average deg coeff for tests A-D, interpolated between design flow temps
    ///
    /// Arguments:
    /// * `flow_temp` - flow temp in K
    pub fn average_degradation_coeff(&self, flow_temp: f64) -> f64 {
        if self.dsgn_flow_temps.len() == 1 {
            return self.average_deg_coeff[0];
        }

        let flow_temp = kelvin_to_celsius(flow_temp);
        interp(
            &self
                .dsgn_flow_temps
                .iter()
                .map(|d| d.0)
                .collect::<Vec<f64>>(),
            &self.average_deg_coeff,
            flow_temp,
        )
    }

    /// Return average capacity for tests A-D, interpolated between design flow temps
    ///
    /// Arguments:
    /// * `flow_temp` - flow temp in K
    pub fn average_capacity(&self, flow_temp: f64) -> f64 {
        if self.dsgn_flow_temps.len() == 1 {
            return self.average_cap[0];
        }

        let flow_temp = kelvin_to_celsius(flow_temp);
        interp(
            &self
                .dsgn_flow_temps
                .iter()
                .map(|d| d.0)
                .collect::<Vec<f64>>(),
            &self.average_cap,
            flow_temp,
        )
    }

    /// Return temperature spread under test conditions, interpolated between design flow temps
    ///
    /// Arguments:
    /// * `flow_temp` - flow temp in K
    pub fn temp_spread_test_conditions(&self, flow_temp: f64) -> f64 {
        if self.dsgn_flow_temps.len() == 1 {
            return self.temp_spread_test_conditions[0];
        }

        let flow_temp = kelvin_to_celsius(flow_temp);
        interp(
            &self
                .dsgn_flow_temps
                .iter()
                .map(|d| d.0)
                .collect::<Vec<f64>>(),
            &self.temp_spread_test_conditions,
            flow_temp,
        )
    }

    fn find_test_record_index(&self, test_condition: &str, dsgn_flow_temp: f64) -> Option<usize> {
        if test_condition == "cld" {
            return Some(0);
        }

        self.test_data[&OrderedFloat(dsgn_flow_temp)]
            .iter()
            .position(|test_record| &test_record.test_letter == test_condition)
    }

    /// Return value at specified test condition, interpolated between design flow temps
    fn data_at_test_condition(
        &self,
        data_item_name: DatumItem,
        test_condition: &str,
        flow_temp: f64,
    ) -> f64 {
        if self.dsgn_flow_temps.len() == 1 {
            let idx = self
                .find_test_record_index(test_condition, self.dsgn_flow_temps[0].0)
                .expect("Expected a test condition to be found in the test data for heat pumps");
            return self.test_data[&self.dsgn_flow_temps[0]][idx].data_item(&data_item_name);
        }

        let data_list = &self
            .dsgn_flow_temps
            .iter()
            .map(|dsgn_flow_temp| {
                let idx = self
                    .find_test_record_index(test_condition, dsgn_flow_temp.0)
                    .expect(
                        "Expected a test condition to be found in the test data for heat pumps",
                    );
                self.test_data[dsgn_flow_temp][idx].data_item(&data_item_name)
            })
            .collect::<Vec<_>>();

        let flow_temp = kelvin_to_celsius(flow_temp);
        interp(
            &self
                .dsgn_flow_temps
                .iter()
                .map(|d| d.0)
                .collect::<Vec<f64>>(),
            data_list,
            flow_temp,
        )
    }

    /// Return Carnot CoP at specified test condition (A, B, C, D, F or cld),
    /// interpolated between design flow temps
    fn carnot_cop_at_test_condition(&self, test_condition: &str, flow_temp: f64) -> f64 {
        self.data_at_test_condition(DatumItem::CarnotCop, test_condition, flow_temp)
    }

    /// Return outlet temp, in Kelvin, at specified test condition (A, B, C, D,
    /// F or cld), interpolated between design flow temps.
    fn outlet_temp_at_test_condition(&self, test_condition: &str, flow_temp: f64) -> f64 {
        celsius_to_kelvin(self.data_at_test_condition(
            DatumItem::TempOutlet,
            test_condition,
            flow_temp,
        ))
    }

    /// Return source temp, in Kelvin, at specified test condition (A, B, C, D,
    /// F or cld), interpolated between design flow temps.
    fn source_temp_at_test_condition(&self, test_condition: &str, flow_temp: f64) -> f64 {
        celsius_to_kelvin(self.data_at_test_condition(
            DatumItem::TempSource,
            test_condition,
            flow_temp,
        ))
    }

    /// Return capacity, in kW, at specified test condition (A, B, C, D, F or
    /// cld), interpolated between design flow temps.
    fn capacity_at_test_condition(&self, test_condition: &str, flow_temp: f64) -> f64 {
        self.data_at_test_condition(DatumItem::Capacity, test_condition, flow_temp)
    }

    /// Return load ratio at operating conditions
    fn load_ratio_at_operating_conditions(
        &self,
        flow_temp: f64,
        temp_source: f64,
        carnot_cop_op_cond: f64,
    ) -> f64 {
        let lr_op_cond_list = self
            .dsgn_flow_temps
            .iter()
            .map(|dsgn_flow_temp| {
                let dsgn_flow_temp = celsius_to_kelvin(dsgn_flow_temp.0);
                let temp_output_cld = self.outlet_temp_at_test_condition("cld", dsgn_flow_temp);
                let temp_source_cld = self.source_temp_at_test_condition("cld", dsgn_flow_temp);
                let carnot_cop_cld = self.carnot_cop_at_test_condition("cld", dsgn_flow_temp);

                let lr_op_cond = (carnot_cop_op_cond / carnot_cop_cld)
                    * (temp_output_cld * temp_source / (flow_temp * temp_source_cld)).powf(N_EXER);
                *[1.0, lr_op_cond]
                    .iter()
                    .max_by(|a, b| a.total_cmp(b))
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let flow_temp = kelvin_to_celsius(flow_temp);
        interp(
            &self
                .dsgn_flow_temps
                .iter()
                .map(|temp| temp.0)
                .collect::<Vec<_>>(),
            &lr_op_cond_list,
            flow_temp,
        )
    }

    /// Return test results either side of operating conditions.
    ///
    /// This function returns 6 results:
    /// - Exergy load ratio below operating conditions
    /// - Exergy load ratio above operating conditions
    /// - Exergy efficiency below operating conditions
    /// - Exergy efficiency above operating conditions
    /// - Degradation coeff below operating conditions
    /// - Degradation coeff above operating conditions
    ///
    /// Arguments:
    /// * `flow_temp` - flow temperature, in Kelvin
    /// * `exergy_lr_op_cond` - exergy load ratio at operating conditions
    fn lr_eff_degcoeff_either_side_of_op_cond(
        &self,
        flow_temp: f64,
        exergy_lr_op_cond: f64,
    ) -> (f64, f64, f64, f64, f64, f64) {
        let mut load_ratios_below: Vec<f64> = Default::default();
        let mut load_ratios_above: Vec<f64> = Default::default();
        let mut efficiencies_below: Vec<f64> = Default::default();
        let mut efficiencies_above: Vec<f64> = Default::default();
        let mut degradation_coeffs_below: Vec<f64> = Default::default();
        let mut degradation_coeffs_above: Vec<f64> = Default::default();

        // For each design flow temperature, find load ratios in test data
        // either side of load ratio calculated for operating conditions.
        // Note: Loop over sorted list of design flow temps and then index into
        //     self.__testdata, rather than looping over self.__testdata,
        //     which is unsorted and therefore may populate the lists in the
        //     wrong order.
        for dsgn_flow_temp in &self.dsgn_flow_temps {
            let dsgn_flow_temp_data = &self.test_data[dsgn_flow_temp];
            // Find the first load ratio in the test data that is greater than
            // or equal to than the load ratio at operating conditions - this
            // and the previous load ratio are the values either side of
            // operating conditions.
            let idx = dsgn_flow_temp_data
                .iter()
                .position(|test_record| test_record.theoretical_load_ratio > exergy_lr_op_cond)
                .unwrap_or(dsgn_flow_temp_data.len() - 1);
            // NB. This assertion reflects an assertion in the source Python.
            assert!(idx > 0);

            // Look up correct load ratio and efficiency based on the idx found above
            load_ratios_below.push(dsgn_flow_temp_data[idx - 1].theoretical_load_ratio);
            load_ratios_above.push(dsgn_flow_temp_data[idx].theoretical_load_ratio);
            efficiencies_below.push(dsgn_flow_temp_data[idx - 1].exergetic_eff);
            efficiencies_above.push(dsgn_flow_temp_data[idx].exergetic_eff);
            degradation_coeffs_below.push(dsgn_flow_temp_data[idx - 1].degradation_coefficient);
            degradation_coeffs_above.push(dsgn_flow_temp_data[idx].degradation_coefficient);
        }

        if self.dsgn_flow_temps.len() == 1 {
            return (
                load_ratios_below[0],
                load_ratios_above[0],
                efficiencies_below[0],
                efficiencies_above[0],
                degradation_coeffs_below[0],
                degradation_coeffs_above[0],
            );
        }

        let flow_temp = kelvin_to_celsius(flow_temp);
        let dsgn_temps_for_interp = &self
            .dsgn_flow_temps
            .iter()
            .map(|temp| temp.0)
            .collect::<Vec<_>>();
        let lr_below = interp(dsgn_temps_for_interp, &load_ratios_below, flow_temp);
        let lr_above = interp(dsgn_temps_for_interp, &load_ratios_above, flow_temp);
        let eff_below = interp(dsgn_temps_for_interp, &efficiencies_below, flow_temp);
        let eff_above = interp(dsgn_temps_for_interp, &efficiencies_above, flow_temp);
        let deg_below = interp(dsgn_temps_for_interp, &degradation_coeffs_below, flow_temp);
        let deg_above = interp(dsgn_temps_for_interp, &degradation_coeffs_above, flow_temp);

        (
            lr_below, lr_above, eff_below, eff_above, deg_below, deg_above,
        )
    }

    /// Calculate CoP at operating conditions when heat pump is not air-source
    ///
    /// Arguments:
    /// * `temp_diff_limit_low` - minimum temperature difference between source and sink
    /// * `temp_ext` - external temperature, in Kelvin
    /// * `temp_source` - source temperature, in Kelvin
    /// * `temp_output` - output temperature, in Kelvin
    fn cop_op_cond_if_not_air_source(
        &self,
        temp_diff_limit_low: f64,
        temp_ext: f64,
        temp_source: f64,
        temp_output: f64,
    ) -> f64 {
        // Need to use Celsius here because regression coeffs were calculated
        // using temperature in Celsius
        let temp_ext = kelvin_to_celsius(temp_ext);

        // For each design flow temperature, calculate CoP at operating conditions
        // Note: Loop over sorted list of design flow temps and then index into
        //        self.test_data, rather than looping over self.test_data,
        //        which is unsorted and therefore may populate the lists in the
        //        wrong order.
        let cop_op_cond = self
            .dsgn_flow_temps
            .iter()
            .map(|dsgn_flow_temp| {
                let dsgn_flow_temp_data = &self.test_data[dsgn_flow_temp];
                // Get the source and outlet temperatures from the coldest test record
                let temp_outlet_cld = celsius_to_kelvin(dsgn_flow_temp_data[0].temp_outlet);
                let temp_source_cld = celsius_to_kelvin(dsgn_flow_temp_data[0].temp_source);

                (self.regression_coeffs[dsgn_flow_temp][0]
                    + self.regression_coeffs[dsgn_flow_temp][1] * temp_ext
                    + self.regression_coeffs[dsgn_flow_temp][2] * temp_ext.powi(2))
                    * temp_output
                    * (temp_outlet_cld - temp_source_cld)
                    / (temp_outlet_cld
                        * *[temp_output - temp_source, temp_diff_limit_low]
                            .iter()
                            .max_by(|a, b| a.total_cmp(b))
                            .unwrap())
            })
            .collect::<Vec<_>>();

        if self.dsgn_flow_temps.len() == 1 {
            return cop_op_cond[0];
        }

        // Interpolate between the values found for the different design flow temperatures
        let flow_temp = kelvin_to_celsius(temp_output);
        interp(
            &self
                .dsgn_flow_temps
                .iter()
                .map(|temp| temp.0)
                .collect::<Vec<_>>(),
            &cop_op_cond,
            flow_temp,
        )
    }

    /// Calculate thermal capacity at operating conditions when heat pump is not air-source
    ///        
    /// Arguments:
    /// * `temp_source` - source temperature, in Kelvin
    /// * `temp_output` - output temperature, in Kelvin
    /// * `mod_ctrl` - boolean specifying whether or not the heat has controls
    ///                    capable of varying the output (as opposed to just on/off
    ///                    control)
    fn capacity_op_cond_if_not_air_source(
        &self,
        temp_output: f64,
        temp_source: f64,
        mod_ctrl: bool,
    ) -> f64 {
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
            self.dsgn_flow_temps.iter().map(|dsgn_flow_temp| {
                let dsgn_flow_temp_data = &self.test_data[dsgn_flow_temp];
                // Get the source and outlet temperatures from the coldest test record
                let temp_outlet_cld = celsius_to_kelvin(dsgn_flow_temp_data[0].temp_outlet);
                let temp_source_cld = celsius_to_kelvin(dsgn_flow_temp_data[0].temp_source);
                // Get the thermal capacity from the coldest test record
                let thermal_capacity_cld = dsgn_flow_temp_data[0].capacity;

                if mod_ctrl {
                    thermal_capacity_cld * ((temp_outlet_cld * temp_source) / (temp_output * temp_source_cld)).powf(N_EXER)
                } else {
                    let d_idx = self.find_test_record_index("D", dsgn_flow_temp.0).expect("Expected to find a test record with the condition 'D' within heat pump test data");
                    // Get the source and outlet temperatures for test condition D
                    let temp_outlet_d = celsius_to_kelvin(dsgn_flow_temp_data[d_idx].temp_outlet);
                    let temp_source_d = celsius_to_kelvin(dsgn_flow_temp_data[d_idx].temp_source);
                    // Get the thermal capacity for test condition D
                    let thermal_capacity_d = dsgn_flow_temp_data[d_idx].capacity;

                    let temp_diff_cld = temp_outlet_cld - temp_source_cld;
                    let temp_diff_d = temp_outlet_d - temp_source_d;
                    let temp_diff_op_cond = temp_output - temp_source;

                    thermal_capacity_cld + (thermal_capacity_d - thermal_capacity_cld) * ((temp_diff_cld - temp_diff_op_cond) / (temp_diff_cld - temp_diff_d))
                }
            }).collect::<Vec<_>>();

        // Interpolate between the values found for the different design flow temperatures
        let flow_temp = kelvin_to_celsius(temp_output);
        interp(
            &self
                .dsgn_flow_temps
                .iter()
                .map(|temp| temp.0)
                .collect::<Vec<_>>(),
            &therm_cap_op_cond,
            flow_temp,
        )
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
    fn temp_spread_correction(
        &self,
        temp_source: f64,
        temp_output: f64,
        temp_diff_evaporator: f64,
        temp_diff_condenser: f64,
        temp_spread_emitter: f64,
    ) -> f64 {
        let temp_spread_correction_list = self
            .dsgn_flow_temps
            .iter()
            .enumerate()
            .map(|(i, dsgn_flow_temp)| {
                let temp_spread_test_cond = self.temp_spread_test_conditions[i];
                1. - ((temp_spread_test_cond - temp_spread_emitter) / 2.)
                    / (temp_output - temp_spread_test_cond / 2. + temp_diff_condenser - temp_source
                        + temp_diff_evaporator)
            })
            .collect::<Vec<_>>();

        let flow_temp = kelvin_to_celsius(temp_output);
        interp(
            &self
                .dsgn_flow_temps
                .iter()
                .map(|d| d.0)
                .collect::<Vec<f64>>(),
            &temp_spread_correction_list,
            flow_temp,
        )
    }
}

/// The list average_deg_coeff will be in the same order as the
/// corresponding elements in self.__dsgn_flow_temps. This behaviour
/// is relied upon elsewhere.
fn ave_degradation_coeff(
    dsgn_flow_temps: &[OrderedFloat<f64>],
    test_data: &HashMap<OrderedFloat<f64>, Vec<HeatPumpTestDatum>>,
) -> Vec<f64> {
    dsgn_flow_temps
        .iter()
        .map(|dsgn_flow_temp| {
            test_data[dsgn_flow_temp]
                .iter()
                .filter(|datum| {
                    TEST_LETTERS_NON_BIVALENT.contains(&datum.test_letter.chars().next().unwrap())
                }) // assumption is made here that the first letter of the "test_letter" attribute can be considered the test letter
                .map(|datum| datum.degradation_coefficient)
                .sum::<f64>()
                / TEST_LETTERS_NON_BIVALENT.len() as f64
        })
        .collect()
}

/// The list average_cap will be in the same order as the
/// corresponding elements in self.__dsgn_flow_temps. This behaviour
/// is relied upon elsewhere.
fn ave_capacity(
    dsgn_flow_temps: &[OrderedFloat<f64>],
    test_data: &HashMap<OrderedFloat<f64>, Vec<HeatPumpTestDatum>>,
) -> Vec<f64> {
    dsgn_flow_temps
        .iter()
        .map(|dsgn_flow_temp| {
            test_data[dsgn_flow_temp]
                .iter()
                .filter(|datum| {
                    TEST_LETTERS_NON_BIVALENT.contains(&datum.test_letter.chars().next().unwrap())
                }) // assumption is made here that the first letter of the "test_letter" attribute can be considered the test letter
                .map(|datum| datum.capacity)
                .sum::<f64>()
                / TEST_LETTERS_NON_BIVALENT.len() as f64
        })
        .collect()
}

const TEMP_SPREAD_AT_TEST_CONDITIONS: [(f64, f64); 4] =
    [(20., 5.0), (35., 5.0), (55., 8.0), (65., 10.0)];

/// List temp spread at test conditions for the design flow temps in the test data
fn init_temp_spread_test_conditions(
    dsgn_flow_temps: &[OrderedFloat<f64>],
) -> Result<Vec<f64>, String> {
    dsgn_flow_temps
        .iter()
        .map(|dsgn_flow_temp| {
            TEMP_SPREAD_AT_TEST_CONDITIONS
                .iter()
                .find_map(|(temp_flow, temp_spread)| {
                    if dsgn_flow_temp == temp_flow {
                        Some(*temp_spread)
                    } else {
                        None
                    }
                })
                .ok_or("Unexpected design flow temperature encountered".to_string())
        })
        .collect()
}

fn init_regression_coeffs(
    dsgn_flow_temps: &Vec<OrderedFloat<f64>>,
    test_data: &HashMap<OrderedFloat<f64>, Vec<HeatPumpTestDatum>>,
) -> Result<HashMap<OrderedFloat<f64>, Vec<f64>>, String> {
    let mut regression_coeffs: HashMap<OrderedFloat<f64>, Vec<f64>> = Default::default();
    for dsgn_flow_temp in dsgn_flow_temps {
        let temp_test_list: Vec<f64> = test_data[dsgn_flow_temp]
            .iter()
            .map(|datum| datum.temp_test)
            .collect();
        let cop_list: Vec<f64> = test_data[dsgn_flow_temp]
            .iter()
            .map(|datum| datum.cop)
            .collect();
        regression_coeffs.insert(*dsgn_flow_temp, polyfit(&temp_test_list, &cop_list, 2)?);
    }

    Ok(regression_coeffs)
}

impl PartialEq for HeatPumpTestDatum {
    fn eq(&self, other: &Self) -> bool {
        self.temp_test == other.temp_test && self.design_flow_temp == other.design_flow_temp
    }
}

const TIME_CONSTANT_WATER: f64 = 1560.;

/// An object to represent a water heating service provided by a heat pump to e.g. a cylinder.
///
/// This object contains the parts of the heat pump calculation that are
/// specific to providing hot water.
#[derive(Clone)]
pub struct HeatPumpServiceWater {
    heat_pump: Arc<Mutex<HeatPump>>,
    service_name: String,
    control: Option<Arc<Control>>,
    temp_hot_water_in_k: f64,
    temp_return_feed_in_k: f64,
    temp_limit_upper_in_k: f64,
    cold_feed: Arc<ColdWaterSource>,
}

impl HeatPumpServiceWater {
    pub fn new(
        heat_pump: Arc<Mutex<HeatPump>>,
        service_name: String,
        temp_hot_water_in_c: f64,
        temp_return_feed_in_c: f64,
        temp_limit_upper_in_c: f64,
        cold_feed: Arc<ColdWaterSource>,
        control: Option<Arc<Control>>,
    ) -> Self {
        Self {
            heat_pump,
            service_name,
            control,
            temp_hot_water_in_k: celsius_to_kelvin(temp_hot_water_in_c),
            temp_return_feed_in_k: celsius_to_kelvin(temp_return_feed_in_c),
            temp_limit_upper_in_k: celsius_to_kelvin(temp_limit_upper_in_c),
            cold_feed,
        }
    }

    pub fn is_on(&self, timestep_idx: usize) -> bool {
        match &self.control {
            Some(ctrl) => ctrl.is_on(timestep_idx),
            None => true,
        }
    }

    /// Calculate the maximum energy output of the HP, accounting for time
    /// spent on higher-priority services
    pub fn energy_output_max(
        &mut self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simulation_time_iteration.index) {
            return 0.0;
        }

        self.heat_pump.lock().unwrap().energy_output_max(
            self.temp_hot_water_in_k,
            self.temp_return_feed_in_k,
            simulation_time_iteration,
        )
    }

    /// Demand energy (in kWh) from the heat pump
    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        let temp_cold_water =
            celsius_to_kelvin(self.cold_feed.temperature(simulation_time_iteration.index));

        let service_on = self.is_on(simulation_time_iteration.index);

        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        self.heat_pump.lock().unwrap().demand_energy(
            &self.service_name,
            &ServiceType::Water,
            energy_demand,
            self.temp_hot_water_in_k,
            self.temp_return_feed_in_k,
            self.temp_limit_upper_in_k,
            TIME_CONSTANT_WATER,
            service_on,
            simulation_time_iteration,
            None,
            Some(temp_cold_water),
        )
    }
}

const TIME_CONSTANT_SPACE_WATER: f64 = 1370.;
const TIME_CONSTANT_SPACE_AIR: f64 = 120.;

#[derive(Clone)]
pub struct HeatPumpServiceSpace {
    heat_pump: Arc<Mutex<HeatPump>>,
    service_name: String,
    control: Arc<Control>,
    temp_limit_upper_in_k: f64,
    temp_diff_emit_dsgn: f64,
}

/// An object to represent a space heating service provided by a heat pump to e.g. radiators.
///
/// This object contains the parts of the heat pump calculation that are
/// specific to providing space heating.
impl HeatPumpServiceSpace {
    /// Arguments:
    /// * `heat_pump` - reference to the HeatPump object providing the service
    /// * `service_name` - name of the service demanding energy from the heat pump
    /// * `temp_limit_upper` - upper operating limit for temperature, in deg C
    /// * `temp_diff_emit_dsgn` - design temperature difference across the emitters, in deg C or K
    /// * `control` - reference to a control object which must implement is_on() and setpnt() funcs
    pub fn new(
        heat_pump: Arc<Mutex<HeatPump>>,
        service_name: String,
        temp_limit_upper_in_c: f64,
        temp_diff_emit_dsgn: f64,
        control: Arc<Control>,
    ) -> Self {
        Self {
            heat_pump,
            service_name,
            control,
            temp_limit_upper_in_k: celsius_to_kelvin(temp_limit_upper_in_c),
            temp_diff_emit_dsgn,
        }
    }

    pub fn is_on(&self, timestep_idx: usize) -> bool {
        self.control.is_on(timestep_idx)
    }

    pub fn temp_setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        per_control!(&self.control.as_ref(), ctrl => { ctrl.setpnt(simulation_time_iteration) })
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        per_control!(&self.control.as_ref(), ctrl => { <_ as ControlBehaviour>::in_required_period(ctrl, simulation_time_iteration) })
    }

    /// Calculate the maximum energy output of the HP, accounting for time
    /// spent on higher-priority services
    pub fn energy_output_max(
        &self,
        temp_output: f64,
        temp_return_feed: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simulation_time_iteration.index) {
            return 0.0;
        }

        let temp_output = celsius_to_kelvin(temp_output);
        self.heat_pump.lock().unwrap().energy_output_max(
            temp_output,
            temp_return_feed,
            simulation_time_iteration,
        )
    }

    /// Demand energy (in kWh) from the heat pump
    ///
    /// Arguments:
    /// * `energy_demand` - space heating energy demand, in kWh
    /// * `temp_flow` - flow temperature for emitters, in deg C
    /// * `temp_return` - return temperature for emitters, in deg C
    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Result<f64, PoisonError<MutexGuard<HeatPump>>> {
        let service_on = self.is_on(simulation_time_iteration.index);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        self.heat_pump.lock().map(|mut pump| {
            let time_constant_for_service = match pump.sink_type {
                HeatPumpSinkType::Water => TIME_CONSTANT_SPACE_WATER,
                HeatPumpSinkType::Air => TIME_CONSTANT_SPACE_AIR,
            };
            pump.demand_energy(
                &self.service_name,
                &ServiceType::Space,
                energy_demand,
                celsius_to_kelvin(temp_flow),
                celsius_to_kelvin(temp_return),
                self.temp_limit_upper_in_k,
                time_constant_for_service,
                service_on,
                simulation_time_iteration,
                Some(TempSpreadCorrectionArg::Callable(
                    self.temp_spread_correction_fn(),
                )),
                None,
            )
        })
    }

    /// Return the cumulative running time and throughput factor (exhaust air HPs only)
    pub fn running_time_throughput_factor(
        &self,
        space_heat_running_time_cumulative: f64,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> (f64, f64) {
        let service_on = self.is_on(simulation_time_iteration.index);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        let time_constant_for_service = match self.heat_pump.lock().unwrap().sink_type {
            HeatPumpSinkType::Water => TIME_CONSTANT_SPACE_WATER,
            HeatPumpSinkType::Air => TIME_CONSTANT_SPACE_AIR,
        };

        self.heat_pump
            .lock()
            .unwrap()
            .running_time_throughput_factor(
                space_heat_running_time_cumulative,
                &self.service_name,
                &ServiceType::Space,
                energy_demand,
                celsius_to_kelvin(temp_flow),
                celsius_to_kelvin(temp_return),
                self.temp_limit_upper_in_k,
                time_constant_for_service,
                service_on,
                simulation_time_iteration,
                Some(TempSpreadCorrectionArg::Callable(
                    self.temp_spread_correction_fn(),
                )),
            )
    }

    fn temp_spread_correction_fn(&self) -> Box<dyn FnOnce(f64, f64) -> f64> {
        // Average temperature difference between heat transfer medium and
        // refrigerant in condenser
        let temp_diff_condenser = 5.0;

        // Average temperature difference between heat transfer medium and
        // refrigerant in evaporator
        let temp_diff_evaporator = match &self.heat_pump.lock().unwrap().source_type {
            t if t.source_fluid_is_air() => 15.0,
            t if t.source_fluid_is_water() => 10.0,
            _ => panic!("impossible heat pump source type encountered"),
        };

        let test_data = self.heat_pump.lock().unwrap().test_data.clone();
        let temp_diff_emit_dsgn = self.temp_diff_emit_dsgn;

        Box::new(move |temp_output, temp_source| {
            test_data.temp_spread_correction(
                temp_source,
                temp_output,
                temp_diff_evaporator,
                temp_diff_condenser,
                temp_diff_emit_dsgn,
            )
        })
    }
}

/// An object to represent a warm air space heating service provided by a heat pump.
///
///    This object contains the parts of the heat pump calculation that are
///    specific to providing space heating via warm air.
#[derive(Clone)]
pub struct HeatPumpServiceSpaceWarmAir {
    heat_pump: Arc<Mutex<HeatPump>>,
    service_name: String,
    control: Arc<Control>,
    temp_limit_upper_in_k: f64,
    temp_diff_emit_dsgn: f64,
    frac_convective: f64,
    temp_flow: f64,
    temp_return: f64,
}

impl HeatPumpServiceSpaceWarmAir {
    pub fn new(
        heat_pump: Arc<Mutex<HeatPump>>,
        service_name: String,
        temp_diff_emit_dsgn: f64,
        control: Arc<Control>,
        temp_flow: f64,
        frac_convective: f64,
    ) -> Self {
        let temp_limit_upper_in_c = temp_flow;

        Self {
            heat_pump,
            service_name,
            control,
            temp_limit_upper_in_k: celsius_to_kelvin(temp_limit_upper_in_c),
            temp_diff_emit_dsgn,
            frac_convective,
            temp_flow,
            temp_return: temp_flow,
        }
    }

    pub fn is_on(&self, timestep_idx: usize) -> bool {
        self.control.is_on(timestep_idx)
    }

    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Result<f64, PoisonError<MutexGuard<HeatPump>>> {
        let temp_flow = self.temp_flow;
        let temp_return = self.temp_return;

        let service_on = self.is_on(simulation_time_iteration.index);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        self.heat_pump.lock().map(|mut pump| {
            let time_constant_for_service = match pump.sink_type {
                HeatPumpSinkType::Water => TIME_CONSTANT_SPACE_WATER,
                HeatPumpSinkType::Air => TIME_CONSTANT_SPACE_AIR,
            };
            pump.demand_energy(
                &self.service_name,
                &ServiceType::Space,
                energy_demand,
                celsius_to_kelvin(temp_flow),
                celsius_to_kelvin(temp_return),
                self.temp_limit_upper_in_k,
                time_constant_for_service,
                service_on,
                simulation_time_iteration,
                Some(TempSpreadCorrectionArg::Callable(
                    self.temp_spread_correction_fn(),
                )),
                None,
            )
        })
    }

    pub fn running_time_throughput_factor(
        &self,
        energy_demand: f64,
        space_heat_running_time_cumulative: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> (f64, f64) {
        let temp_flow = self.temp_flow;
        let temp_return = self.temp_return;

        let service_on = self.is_on(simulation_time_iteration.index);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        let time_constant_for_service = match self.heat_pump.lock().unwrap().sink_type {
            HeatPumpSinkType::Water => TIME_CONSTANT_SPACE_WATER,
            HeatPumpSinkType::Air => TIME_CONSTANT_SPACE_AIR,
        };

        self.heat_pump
            .lock()
            .unwrap()
            .running_time_throughput_factor(
                space_heat_running_time_cumulative,
                &self.service_name,
                &ServiceType::Space,
                energy_demand,
                celsius_to_kelvin(temp_flow),
                celsius_to_kelvin(temp_return),
                self.temp_limit_upper_in_k,
                time_constant_for_service,
                service_on,
                simulation_time_iteration,
                Some(TempSpreadCorrectionArg::Callable(
                    self.temp_spread_correction_fn(),
                )),
            )
    }

    /// yes this is copy-pasted from the SpaceService
    fn temp_spread_correction_fn(&self) -> Box<dyn FnOnce(f64, f64) -> f64> {
        // Average temperature difference between heat transfer medium and
        // refrigerant in condenser
        let temp_diff_condenser = 5.0;

        // Average temperature difference between heat transfer medium and
        // refrigerant in evaporator
        let temp_diff_evaporator = match &self.heat_pump.lock().unwrap().source_type {
            t if t.source_fluid_is_air() => 15.0,
            t if t.source_fluid_is_water() => 10.0,
            _ => panic!("impossible heat pump source type encountered"),
        };

        let test_data = self.heat_pump.lock().unwrap().test_data.clone();
        let temp_diff_emit_dsgn = self.temp_diff_emit_dsgn;

        Box::new(move |temp_output, temp_source| {
            test_data.temp_spread_correction(
                temp_source,
                temp_output,
                temp_diff_evaporator,
                temp_diff_condenser,
                temp_diff_emit_dsgn,
            )
        })
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
#[derive(Clone)]
pub struct HeatPump {
    // energy supply
    simulation_timestep: f64,
    external_conditions: Arc<ExternalConditions>,
    throughput_exhaust_air: f64,
    // litres/ second
    service_results: Arc<Mutex<Vec<ServiceResult>>>,
    total_time_running_current_timestep: f64,
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
    power_source_circ_pump: f64,
    power_standby: f64,
    power_crankcase_heater_mode: f64,
    power_off_mode: f64,
    power_max_backup: f64,
    temp_distribution_heat_network: Option<f64>,
    // energy supply for heat networks
    overvent_ratio: f64,
    test_data: HeatPumpTestData,
    temp_min_modulation_rate_low: Option<f64>,
    min_modulation_rate_low: Option<f64>,
    temp_min_modulation_rate_high: Option<f64>,
    min_modulation_rate_55: Option<f64>,
    detailed_results: Option<Vec<Arc<Mutex<Vec<ServiceResult>>>>>,
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
    /// * `throughput_exhaust_air` - throughput (litres / second) of exhaust air
    /// * `heat_network` - reference to EnergySupply object representing heat network
    ///                        (for HPs that use heat network as heat source)
    /// * `output_detailed_results` - if true, save detailed results from each timestep
    ///                                   for later reporting
    ///
    /// Other variables:
    ///        energy_supply_connections
    ///            -- dictionary with service name strings as keys and corresponding
    ///               EnergySupplyConnection objects as values
    ///        energy_supply_connection_aux -- EnergySupplyConnection object for auxiliary energy
    ///        test_data -- HeatPumpTestData object
    pub fn new(
        heat_pump_input: &HeatSourceWetDetails,
        simulation_timestep: f64,
        external_conditions: Arc<ExternalConditions>,
        throughput_exhaust_air: Option<f64>,
        // heat_network: Option<Arc<HeatNetwork>>,
        output_detailed_results: bool,
    ) -> Result<Self, String> {
        let service_results = Default::default();
        let total_time_running_current_timestep = Default::default();
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
                power_source_circ_pump,
                power_standby,
                power_crankcase_heater,
                power_off,
                power_max_backup,
                temp_distribution_heat_network,
                test_data,
                min_modulation_rate_20,
                min_modulation_rate_35,
                min_modulation_rate_55,
                ..
            } => (
                *source_type,
                *sink_type,
                *backup_control_type,
                *time_delay_backup,
                *modulating_control,
                *time_constant_onoff_operation,
                *temp_return_feed_max,
                celsius_to_kelvin(*temp_lower_operating_limit),
                *min_temp_diff_flow_return_for_hp_to_operate,
                *var_flow_temp_ctrl_during_test,
                *power_heating_circ_pump,
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
            ),
            _ => panic!(
                "Expected heat source details passed to heat pump to be of heat pump variant"
            ),
        };
        let time_delay_backup = match backup_ctrl {
            HeatPumpBackupControlType::None => None,
            _ => Some(time_delay_backup),
        };
        let temp_return_feed_max = match sink_type {
            HeatPumpSinkType::Air => None,
            _ => Some(celsius_to_kelvin(temp_return_feed_max)),
        };
        let temp_distribution_heat_network = match source_type {
            HeatPumpSourceType::HeatNetwork => *temp_distribution_heat_network,
            _ => None,
        };

        // Exhaust air HP requires different/additional initialisation, which is implemented here
        let overvent_ratio = if source_type.is_exhaust_air() {
            let (lowest_air_flow_rate_in_test_data, _) =
                interpolate_exhaust_air_heat_pump_test_data(
                    throughput_exhaust_air
                        .expect("expected throughput_exhaust_air to have been set here"),
                    &test_data,
                )?;
            *[
                1.0,
                lowest_air_flow_rate_in_test_data / throughput_exhaust_air.unwrap(),
            ]
            .iter()
            .max_by(|a, b| a.total_cmp(b))
            .unwrap()
        } else {
            1.0
        };

        // in the Python there is a check here about an air flow rate being present - decided
        // any validation should really happen upstream of here

        let test_data = HeatPumpTestData::new(test_data)?;

        let (
            temp_min_modulation_rate_low,
            min_modulation_rate_low,
            temp_min_modulation_rate_high,
            min_modulation_rate_55,
        ) = if modulating_ctrl {
            let (temp_min_modulation_rate_low, min_modulation_rate_low) = match sink_type {
                HeatPumpSinkType::Air => (
                    Some(20.),
                    Some(
                        min_modulation_rate_20
                            .ok_or("expected a min_modulation_rate_20 provided".to_string())?,
                    ),
                ),
                HeatPumpSinkType::Water => (
                    Some(35.),
                    Some(
                        min_modulation_rate_35
                            .ok_or("expected a min_modulation_rate_35 provided".to_string())?,
                    ),
                ),
            };
            let (temp_min_modulation_rate_high, min_modulation_rate_55) =
                if test_data.dsgn_flow_temps.contains(&OrderedFloat(55.)) {
                    (
                        Some(55.),
                        Some(
                            min_modulation_rate_55
                                .ok_or("expected a min_modulation_rate_55 provided".to_string())?,
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

        Ok(Self {
            simulation_timestep,
            external_conditions,
            throughput_exhaust_air: throughput_exhaust_air.unwrap(), // this may not be optional
            service_results,
            total_time_running_current_timestep,
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
            power_heating_circ_pump,
            power_source_circ_pump,
            power_standby,
            power_crankcase_heater_mode,
            power_off_mode,
            power_max_backup,
            temp_distribution_heat_network,
            overvent_ratio,
            test_data,
            temp_min_modulation_rate_low,
            min_modulation_rate_low,
            temp_min_modulation_rate_high,
            min_modulation_rate_55,
            detailed_results,
        })
    }

    pub fn source_is_exhaust_air(&self) -> bool {
        self.source_type.is_exhaust_air()
    }

    // TODO implement create_service_connection method

    // TODO implement create service methods

    pub fn create_service_hot_water(
        heat_pump: Arc<Mutex<Self>>,
        service_name: String,
        temp_hot_water_in_c: f64,
        temp_return_feed_in_c: f64,
        temp_limit_upper_in_c: f64,
        cold_feed: Arc<ColdWaterSource>,
        control: Option<Arc<Control>>,
    ) -> HeatPumpServiceWater {
        HeatPumpServiceWater::new(
            heat_pump,
            service_name,
            temp_hot_water_in_c,
            temp_return_feed_in_c,
            temp_limit_upper_in_c,
            cold_feed,
            control,
        )
    }

    pub fn create_service_space_heating(
        heat_pump: Arc<Mutex<Self>>,
        service_name: String,
        temp_limit_upper_in_c: f64,
        temp_diff_emit_dsgn: f64,
        control: Arc<Control>,
    ) -> HeatPumpServiceSpace {
        HeatPumpServiceSpace::new(
            heat_pump,
            service_name,
            temp_limit_upper_in_c,
            temp_diff_emit_dsgn,
            control,
        )
    }

    pub fn create_service_space_heating_warm_air(
        self,
        service_name: String,
        control: Arc<Control>,
        frac_convective: f64,
    ) -> Result<HeatPumpServiceSpaceWarmAir, &'static str> {
        if !matches!(&self.sink_type, HeatPumpSinkType::Air) {
            return Err("Warm air space heating service requires heat pump with sink type Air");
        }

        // Use low temperature test data for space heating - set flow temp such
        // that it matches the one used in the test
        let temp_flow = self.test_data.dsgn_flow_temps[0].0;

        // Design temperature difference across the emitters, in deg C or K
        let temp_diff_emit_dsgn = max_of_2(
            temp_flow / 7.,
            self.test_data.temp_spread_test_conditions(temp_flow),
        );

        Ok(HeatPumpServiceSpaceWarmAir::new(
            Arc::new(Mutex::new(self)),
            service_name,
            temp_diff_emit_dsgn,
            control,
            temp_flow,
            frac_convective,
        ))
    }

    /// Get source temp according to rules in CALCM-01 - DAHPSE - V2.0_DRAFT13, 3.1.1
    fn get_temp_source(&self, simulation_time_iteration: &SimulationTimeIteration) -> f64 {
        let temp_source = match self.source_type {
            HeatPumpSourceType::Ground => {
                let temp_ext = self.external_conditions.air_temp(simulation_time_iteration);
                *[
                    0.,
                    *[8., temp_ext * 0.25806 + 2.8387]
                        .iter()
                        .max_by(|a, b| a.total_cmp(b).reverse())
                        .unwrap(),
                ]
                    .iter()
                    .max_by(|a, b| a.total_cmp(b))
                    .unwrap()
            }
            HeatPumpSourceType::OutsideAir => self.external_conditions.air_temp(simulation_time_iteration),
            HeatPumpSourceType::ExhaustAirMEV => 20.0,
            HeatPumpSourceType::ExhaustAirMVHR => 20.0,
            HeatPumpSourceType::ExhaustAirMixed => {
                todo!() // TODO is from Python
            }
            HeatPumpSourceType::WaterGround => {
                todo!() // TODO is from Python
            }
            HeatPumpSourceType::WaterSurface => {
                todo!() // TODO is from Python
            }
            HeatPumpSourceType::HeatNetwork => self.temp_distribution_heat_network.expect("expected a temperature to be available for a heat network if this is the source type"),
        };

        celsius_to_kelvin(temp_source)
    }

    /// Calculate the thermal capacity of the heat pump at operating conditions
    ///
    /// Based on CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.4
    fn thermal_capacity_op_cond(&self, temp_output: f64, temp_source: f64) -> f64 {
        if !matches!(self.source_type, HeatPumpSourceType::OutsideAir)
            && !self.var_flow_temp_ctrl_during_test
        {
            self.test_data.average_capacity(temp_output)
        } else {
            self.test_data.capacity_op_cond_if_not_air_source(
                temp_output,
                temp_source,
                self.modulating_ctrl,
            )
        }
    }

    /// Calculate the maximum energy output of the HP, accounting for time
    /// spent on higher-priority services
    ///
    /// Note: Call via a HeatPumpService object, not directly.
    pub fn energy_output_max(
        &self,
        temp_output: f64,
        temp_return_feed: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        let timestep = self.simulation_timestep;
        let time_available = timestep - self.total_time_running_current_timestep;
        let temp_source = self.get_temp_source(simulation_time_iteration);

        let power_max_hp =
            if self.outside_operating_limits(temp_return_feed, simulation_time_iteration) {
                0.0
            } else {
                self.thermal_capacity_op_cond(temp_output, temp_source)
            };

        let power_max = if matches!(self.backup_ctrl, HeatPumpBackupControlType::None)
            || !self.backup_heater_delay_time_elapsed()
        {
            power_max_hp
        } else if matches!(self.backup_ctrl, HeatPumpBackupControlType::TopUp) {
            power_max_hp + self.power_max_backup
        } else {
            // arm for HeatPumpBackupControlType::Substitute as Rust requires assigning
            // conditional to be explicitly exhaustive and we can't do this using a match here
            *[power_max_hp, self.power_max_backup]
                .iter()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap()
        };

        power_max * time_available
    }

    /// Calculate CoP and degradation coefficient at operating conditions
    fn cop_deg_coeff_op_cond(
        &self,
        service_type: &ServiceType,
        temp_output: f64,
        temp_source: f64,
        temp_spread_correction: TempSpreadCorrectionArg,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> (f64, f64) {
        let temp_spread_correction_factor = match temp_spread_correction {
            TempSpreadCorrectionArg::Float(correction) => correction,
            TempSpreadCorrectionArg::Callable(callable) => callable(temp_output, temp_source),
        };

        if !matches!(self.source_type, HeatPumpSourceType::OutsideAir)
            && !self.var_flow_temp_ctrl_during_test
        {
            let cop_op_cond = temp_spread_correction_factor
                * self.test_data.cop_op_cond_if_not_air_source(
                    HEAT_PUMP_TEMP_DIFF_LIMIT_LOW,
                    self.external_conditions.air_temp(simulation_time_iteration), // TODO Python uses .temperature() method here although need to check if this is a bug
                    temp_source,
                    temp_output,
                );
            let deg_coeff_op_cond = self.test_data.average_degradation_coeff(temp_output);

            (cop_op_cond, deg_coeff_op_cond)
        } else {
            let carnot_cop_op_cond = carnot_cop(
                temp_source,
                temp_output,
                Some(HEAT_PUMP_TEMP_DIFF_LIMIT_LOW),
            );
            // Get exergy load ratio at operating conditions and exergy load ratio,
            // exergy efficiency and degradation coeff at test conditions above and
            // below operating conditions
            let lr_op_cond = self.test_data.load_ratio_at_operating_conditions(
                temp_output,
                temp_source,
                carnot_cop_op_cond,
            );
            let (lr_below, lr_above, eff_below, eff_above, deg_coeff_below, deg_coeff_above) = self
                .test_data
                .lr_eff_degcoeff_either_side_of_op_cond(temp_output, lr_op_cond);

            // CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.4
            // Get exergy efficiency by interpolating between figures above and
            // below operating conditions
            let exer_eff_op_cond = eff_below
                + (eff_below - eff_above) * (lr_op_cond - lr_below) / (lr_below - lr_above);

            // CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.5
            // Note: DAHPSE method document section 4.5.5 doesn't have
            // temp_spread_correction_factor in formula below. However, section 4.5.7
            // states that the correction factor is to be applied to the CoP.
            let cop_op_cond = *[
                1.0,
                exer_eff_op_cond * carnot_cop_op_cond * temp_spread_correction_factor,
            ]
            .iter()
            .max_by(|a, b| a.total_cmp(b))
            .unwrap();

            let (limit_upper, limit_lower) = if matches!(self.sink_type, HeatPumpSinkType::Air)
                && !matches!(service_type, ServiceType::Water)
            {
                (0.25, 0.0)
            } else {
                (1.0, 0.9)
            };

            let deg_coeff_op_cond = if lr_below == lr_above {
                deg_coeff_below
            } else {
                deg_coeff_below
                    + (deg_coeff_below - deg_coeff_above) * (lr_op_cond - lr_below)
                        / (lr_below - lr_above)
            };

            let deg_coeff_op_cond = *[
                *[deg_coeff_op_cond, limit_upper]
                    .iter()
                    .max_by(|a, b| a.total_cmp(b).reverse())
                    .unwrap(),
                limit_lower,
            ]
            .iter()
            .max_by(|a, b| a.total_cmp(b))
            .unwrap();

            (cop_op_cond, deg_coeff_op_cond)
        }
    }

    /// Calculate energy output limited by upper temperature
    fn energy_output_limited(
        &self,
        energy_output_required: f64,
        temp_output: f64,
        temp_used_for_scaling: f64,
        temp_limit_upper: f64,
    ) -> f64 {
        if temp_output > temp_limit_upper {
            if temp_output == temp_used_for_scaling {
                energy_output_required
            } else {
                if (temp_limit_upper - temp_used_for_scaling) >= self.temp_diff_flow_return_min {
                    energy_output_required * (temp_limit_upper - temp_used_for_scaling)
                        / (temp_output - temp_used_for_scaling)
                } else {
                    0.
                }
            }
        } else {
            energy_output_required
        }
    }

    /// Check if backup heater is available or still in delay period
    fn backup_heater_delay_time_elapsed(&self) -> bool {
        self.time_running_continuous
            >= self
                .time_delay_backup
                .expect("time delay backup was expected to be set")
    }

    /// Check if heat pump is outside operating limits
    fn outside_operating_limits(
        &self,
        temp_return_feed: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> bool {
        let temp_source = self.get_temp_source(simulation_time_iteration);
        let below_min_ext_temp = temp_source <= self.temp_lower_op_limit;

        let above_temp_return_feed_max = match self.sink_type {
            HeatPumpSinkType::Water => {
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
    ) -> bool {
        let timestep = self.simulation_timestep;

        // For top-up backup heater, use backup if delay time has elapsed.
        // For substitute backup heater, use backup if delay time has elapsed and
        // backup heater can provide more energy than heat pump. This assumption
        // is required to make the maximum energy output of the system
        // predictable before the demand is known.
        if (matches!(self.backup_ctrl, HeatPumpBackupControlType::TopUp)
            && self.backup_heater_delay_time_elapsed())
            || (matches!(self.backup_ctrl, HeatPumpBackupControlType::Substitute)
                && self.backup_heater_delay_time_elapsed()
                && self.power_max_backup > thermal_capacity_op_cond)
        {
            energy_output_required > thermal_capacity_op_cond * timestep
        } else {
            false
        }
    }

    /// Evaluate boolean conditions that may trigger backup heater
    fn use_backup_heater_only(
        &self,
        energy_output_required: f64,
        temp_return_feed: f64,
        thermal_capacity_op_cond: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> bool {
        let outside_operating_limits =
            self.outside_operating_limits(temp_return_feed, simulation_time_iteration);
        let inadequate_capacity =
            self.inadequate_capacity(energy_output_required, thermal_capacity_op_cond);

        !matches!(self.backup_ctrl, HeatPumpBackupControlType::None)
            && (outside_operating_limits
                || (inadequate_capacity
                    && matches!(self.backup_ctrl, HeatPumpBackupControlType::Substitute)))
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
    ///     self.total_time_running_current_timestep
    fn run_demand_energy_calc(
        &self,
        service_name: &str,
        service_type: &ServiceType,
        energy_output_required: f64,
        temp_output: f64,
        temp_return_feed: f64,
        temp_limit_upper: f64,
        time_constant_for_service: f64,
        service_on: bool,
        simulation_time_iteration: &SimulationTimeIteration,
        temp_spread_correction: Option<TempSpreadCorrectionArg>,
        temp_used_for_scaling: Option<f64>,
        additional_time_unavailable: Option<f64>,
    ) -> HeatPumpEnergyCalculation {
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

        let temp_source = self.get_temp_source(simulation_time_iteration);
        // From here onwards, output temp to be used is subject to the upper limit
        let temp_output = min_of_2(temp_output, temp_limit_upper);

        // Get thermal capacity, CoP and degradation coeff at operating conditions
        let thermal_capacity_op_cond = self.thermal_capacity_op_cond(temp_output, temp_source);
        let (cop_op_cond, deg_coeff_op_cond) = self.cop_deg_coeff_op_cond(
            service_type,
            temp_output,
            temp_source,
            temp_spread_correction,
            simulation_time_iteration,
        );

        // Calculate running time of HP
        let time_required = energy_output_limited / thermal_capacity_op_cond;
        let time_available =
            timestep - self.total_time_running_current_timestep - additional_time_unavailable;
        let time_running_current_service = min_of_2(time_required, time_available);

        let temp_min_modulation_rate_low = self
            .temp_min_modulation_rate_low
            .expect("temp_min_modulation_rate_low expected to have been provided");
        let temp_min_modulation_rate_high = self
            .temp_min_modulation_rate_high
            .expect("temp_min_modulation_rate_high expected to have been provided");
        let min_modulation_rate_low = self
            .min_modulation_rate_low
            .expect("A min_modulation_rate_low was expected to have been set");

        // Calculate load ratio
        let load_ratio = time_running_current_service / timestep;
        let load_ratio_continuous_min = if self.modulating_ctrl {
            if self.test_data.dsgn_flow_temps.contains(&OrderedFloat(55.)) {
                (min_of_2(
                    max_of_2(temp_output, temp_min_modulation_rate_low),
                    temp_min_modulation_rate_high,
                ) - temp_min_modulation_rate_low)
                    / (temp_min_modulation_rate_high - temp_min_modulation_rate_low)
                    * min_modulation_rate_low
                    + (1.0
                        - (min_of_2(
                            max_of_2(temp_output, temp_min_modulation_rate_low),
                            temp_min_modulation_rate_high,
                        ) - temp_min_modulation_rate_low))
                        / (temp_min_modulation_rate_high - temp_min_modulation_rate_low)
                        * self
                            .min_modulation_rate_55
                            .expect("A min modulation rate was expected for a 55 air flow temp")
            } else {
                min_modulation_rate_low
            }
        } else {
            // On/off heat pump cannot modulate below maximum power
            1.0
        };

        // Determine whether or not HP is operating in on/off mode
        let hp_operating_in_onoff_mode = load_ratio > 0.0 && load_ratio < load_ratio_continuous_min;

        let compressor_power_full_load = thermal_capacity_op_cond / cop_op_cond;

        // CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.10, step 1:
        let compressor_power_min_load = compressor_power_full_load * load_ratio_continuous_min;
        // CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.5.10, step 2:
        let power_used_due_to_inertia_effects = if load_ratio >= load_ratio_continuous_min {
            0.0
        } else {
            compressor_power_min_load
                * self.time_constant_onoff_operation
                * load_ratio
                * (1. - load_ratio)
                / time_constant_for_service
        };

        let use_backup_heater_only = self.use_backup_heater_only(
            energy_output_required,
            temp_return_feed,
            thermal_capacity_op_cond,
            simulation_time_iteration,
        );

        // Calculate energy delivered by HP and energy input
        let (energy_delivered_hp, energy_input_hp, energy_input_hp_divisor) =
            if use_backup_heater_only || !service_on {
                (0., 0., None)
            } else {
                let energy_delivered_hp = thermal_capacity_op_cond * time_running_current_service;
                let (energy_input_hp, energy_input_hp_divisor) = if hp_operating_in_onoff_mode {
                    let energy_input_hp_divisor = if matches!(service_type, ServiceType::Water)
                        && matches!(self.sink_type, HeatPumpSinkType::Air)
                    {
                        1. - deg_coeff_op_cond * (1. - load_ratio / load_ratio_continuous_min)
                    } else {
                        1.
                    };

                    // Note: Energy_ancillary_when_off should also be included in the
                    // energy input for on/off operation, but at this stage we have
                    // not calculated whether a lower-priority service will run
                    // instead, so this will need to be calculated later and
                    // (energy_ancillary_when_off / eqn_denom) added to the energy
                    // input
                    let energy_input_hp = ((compressor_power_full_load * (1. + HEAT_PUMP_F_AUX)
                        + power_used_due_to_inertia_effects)
                        * time_running_current_service)
                        / energy_input_hp_divisor;

                    (energy_input_hp, Some(energy_input_hp_divisor))
                } else {
                    (energy_delivered_hp / cop_op_cond, None)
                };

                (
                    energy_delivered_hp,
                    energy_input_hp,
                    energy_input_hp_divisor,
                )
            };

        // Calculate energy delivered by backup heater
        let energy_delivered_backup = if matches!(self.backup_ctrl, HeatPumpBackupControlType::None)
            || !self.backup_heater_delay_time_elapsed()
            || !service_on
        {
            0.0
        } else {
            // arm here in Python covered other backup control types explicitly - here we're exhaustively just matching everything else
            max_of_2(
                min_of_2(
                    self.power_max_backup * time_available,
                    energy_output_required - energy_delivered_hp,
                ),
                0.0,
            )
        };

        // Calculate energy input to backup heater
        let energy_input_backup = energy_delivered_backup;

        // Energy used by pumps
        let energy_heating_circ_pump = time_running_current_service * self.power_heating_circ_pump;
        let energy_source_circ_pump = time_running_current_service * self.power_source_circ_pump;

        // Calculate total energy delivered and input
        let energy_delivered_total = energy_delivered_hp + energy_delivered_backup;
        let energy_input_total = energy_input_hp
            + energy_input_backup
            + energy_heating_circ_pump
            + energy_source_circ_pump;

        HeatPumpEnergyCalculation {
            service_name: result_str(service_name),
            service_type: *service_type,
            service_on,
            energy_output_required,
            temp_output,
            temp_source,
            cop_op_cond,
            thermal_capacity_op_cond,
            time_running: time_running_current_service,
            deg_coeff_op_cond,
            compressor_power_min_load,
            load_ratio_continuous_min,
            load_ratio,
            use_backup_heater_only,
            hp_operating_in_onoff_mode,
            energy_input_hp_divisor,
            energy_input_hp,
            energy_delivered_hp,
            energy_input_backup,
            energy_delivered_backup,
            energy_input_total,
            energy_delivered_total,
            energy_heating_circ_pump,
            energy_source_circ_pump,
        }
    }

    /// Calculate energy required by heat pump to satisfy demand for the service indicated.
    ///
    /// Note: Call via a HeatPumpService object, not directly.
    pub fn demand_energy(
        &mut self,
        service_name: &str,
        service_type: &ServiceType,
        energy_output_required: f64,
        temp_output: f64,
        temp_return_feed: f64,
        temp_limit_upper: f64,
        time_constant_for_service: f64,
        service_on: bool,
        simulation_time_iteration: &SimulationTimeIteration,
        temp_spread_correction: Option<TempSpreadCorrectionArg>,
        temp_used_for_scaling: Option<f64>,
    ) -> f64 {
        let service_results = self.run_demand_energy_calc(
            service_name,
            service_type,
            energy_output_required,
            temp_output,
            temp_return_feed,
            temp_limit_upper,
            time_constant_for_service,
            service_on,
            simulation_time_iteration,
            temp_spread_correction,
            temp_used_for_scaling,
            None,
        );

        let time_running = service_results.time_running;
        let energy_delivered_total = service_results.energy_delivered_total;

        // Save results that are needed later (in the timestep_end function)
        if let Ok(mut results) = self.service_results.lock() {
            results.push(ServiceResult::Full(service_results));
        }
        self.total_time_running_current_timestep += time_running;

        // Feed/return results to other modules
        // TODO report energy to energy supply

        energy_delivered_total
    }

    /// Return the cumulative running time and throughput factor (exhaust air HPs only)
    pub fn running_time_throughput_factor(
        &self,
        space_heat_running_time_cumulative: f64,
        service_name: &str,
        service_type: &ServiceType,
        energy_output_required: f64,
        temp_output: f64,
        temp_return_feed: f64,
        temp_limit_upper: f64,
        time_constant_for_service: f64,
        service_on: bool,
        simulation_time_iteration: &SimulationTimeIteration,
        temp_spread_correction: Option<TempSpreadCorrectionArg>,
    ) -> (f64, f64) {
        let timestep = self.simulation_timestep;

        let service_results = self.run_demand_energy_calc(
            service_name,
            service_type,
            energy_output_required,
            temp_output,
            temp_return_feed,
            temp_limit_upper,
            time_constant_for_service,
            service_on,
            simulation_time_iteration,
            temp_spread_correction,
            None,
            Some(space_heat_running_time_cumulative),
        );

        // Add running time for the current space heating service to the cumulative total
        let space_heat_running_time_cumulative =
            space_heat_running_time_cumulative + service_results.time_running;

        // Total running time for calculating throughput factor is time already
        // spent on water heating plus (unsaved) time running space heating. This
        // assumes that water heating service is always calculated before space
        // heating service.
        let total_running_time =
            self.total_time_running_current_timestep + space_heat_running_time_cumulative;

        // Apply overventilation ratio to part of timestep where HP is running
        // to calculate throughput_factor.
        let throughput_factor =
            ((timestep - total_running_time) + self.overvent_ratio * total_running_time) / timestep;

        (total_running_time, throughput_factor)
    }

    fn calc_ancillary_energy(&mut self, timestep: f64, time_remaining_current_timestep: f64) {
        // we need to collect times running separately before the main iteration as we cannot look into the service
        // results while iterating over mutable values
        if let Ok(mut service_results) = self.service_results.lock() {
            let times_running_subsequent_services = service_results
                .iter()
                .filter(|r| matches!(r, ServiceResult::Full(_)))
                .enumerate()
                .map(|(service_no, _)| {
                    service_results[(service_no + 1)..]
                        .iter()
                        .filter(|r| matches!(r, ServiceResult::Full(_)))
                        .map(|result| {
                            if let ServiceResult::Full(result) = result {
                                result.time_running
                            } else {
                                panic!("full calculation expected")
                            }
                        })
                        .sum::<f64>()
                })
                .collect::<Vec<f64>>();

            for (service_no, service_data) in service_results.iter_mut().enumerate() {
                // we always expect full results here, but to be explicit:
                if let ServiceResult::Full(ref mut service_data) = service_data {
                    // Unpack results of previous calculations for this service
                    let _service_name = service_data.service_name.as_str();
                    let service_type = service_data.service_type;
                    let service_on = service_data.service_on;
                    let time_running_current_service = service_data.time_running;
                    let deg_coeff_op_cond = service_data.deg_coeff_op_cond;
                    let compressor_power_min_load = service_data.compressor_power_min_load;
                    let load_ratio_continuous_min = service_data.load_ratio_continuous_min;
                    let load_ratio = service_data.load_ratio;
                    let use_backup_heater_only = service_data.use_backup_heater_only;
                    let hp_operating_in_onoff_mode = service_data.hp_operating_in_onoff_mode;
                    let energy_input_hp_divisor = service_data
                        .energy_input_hp_divisor
                        .expect("expected energy_input_hp_divisor to be set in a test record");

                    let time_running_subsequent_services =
                        times_running_subsequent_services[service_no];

                    #[allow(clippy::neg_cmp_op_on_partial_ord)]
                    let energy_ancillary_when_off = if service_on
                        && time_running_current_service > 0.
                        && !(time_running_subsequent_services > 0.)
                        && !(matches!(self.sink_type, HeatPumpSinkType::Air)
                            && matches!(service_type, ServiceType::Water))
                    {
                        (1. - deg_coeff_op_cond)
                            * (compressor_power_min_load / load_ratio_continuous_min)
                            * max_of_2(
                                (time_remaining_current_timestep
                                    - load_ratio / load_ratio_continuous_min * timestep),
                                0.,
                            )
                    } else {
                        0.
                    };

                    let energy_input_hp = if !use_backup_heater_only && hp_operating_in_onoff_mode {
                        energy_ancillary_when_off / energy_input_hp_divisor
                    } else {
                        0.
                    };

                    // TODO record energy use to energy supply
                    service_data.energy_input_hp += energy_input_hp;
                    service_data.energy_input_total += energy_input_hp;
                }
            }
        }
    }

    /// Calculate auxiliary energy according to CALCM-01 - DAHPSE - V2.0_DRAFT13, section 4.7
    fn calc_auxiliary_energy(
        &self,
        timestep: f64,
        time_remaining_current_timestep: f64,
    ) -> (f64, f64, f64) {
        // Retrieve control settings for this timestep
        let mut heating_profile_on = false;
        let mut water_profile_on = false;
        if let Ok(service_results) = self.service_results.lock() {
            for service_data in service_results.iter() {
                if let ServiceResult::Full(calculation) = service_data {
                    match calculation.service_type {
                        ServiceType::Space => {
                            heating_profile_on = calculation.service_on;
                        }
                        ServiceType::Water => {
                            water_profile_on = calculation.service_on;
                        }
                    }
                }
            }
        }

        // Energy used in standby and crankcase heater mode
        let mut energy_off_mode = 0.0;
        let mut energy_standby = 0.0;
        let mut energy_crankcase_heater_mode = 0.0;
        match (heating_profile_on, water_profile_on) {
            (true, _) => {
                energy_standby = time_remaining_current_timestep * self.power_standby;
                energy_crankcase_heater_mode =
                    time_remaining_current_timestep * self.power_crankcase_heater_mode;
            }
            (false, true) => {
                energy_standby = time_remaining_current_timestep * self.power_standby;
            }
            (false, false) => {
                energy_off_mode = timestep * self.power_off_mode;
            }
        }

        let _energy_aux = energy_standby + energy_crankcase_heater_mode + energy_off_mode;
        // TODO report energy demand to energy supply

        (
            energy_standby,
            energy_crankcase_heater_mode,
            energy_off_mode,
        )
    }

    /// If HP uses heat network as source, calculate energy extracted from heat network
    fn extract_energy_from_source(&self) {
        if let Ok(service_results) = self.service_results.lock() {
            for service_data in service_results.iter() {
                if let ServiceResult::Full(service_data) = service_data {
                    let HeatPumpEnergyCalculation {
                        service_name: _service_name,
                        energy_delivered_hp,
                        energy_input_hp,
                        ..
                    } = service_data;
                    let _energy_extracted_hp = energy_delivered_hp - energy_input_hp;
                    // TODO report energy demand to energy supply connections
                }
            }
        }
    }

    /// Calculations to be done at the end of each timestep
    fn timestep_end(&mut self) {
        let timestep = self.simulation_timestep;
        let time_remaining_current_timestep = timestep - self.total_time_running_current_timestep;

        if time_remaining_current_timestep == 0.0 {
            self.time_running_continuous += self.total_time_running_current_timestep;
        }
        {
            self.time_running_continuous = 0.;
        }

        self.calc_ancillary_energy(timestep, time_remaining_current_timestep);
        let (energy_standby, energy_crankcase_heater_mode, energy_off_mode) =
            self.calc_auxiliary_energy(timestep, time_remaining_current_timestep);

        if matches!(self.source_type, HeatPumpSourceType::HeatNetwork) {
            self.extract_energy_from_source();
        }

        // If detailed results are to be output, save the results from the current timestep
        if let Some(ref mut detailed_results) = self.detailed_results {
            if let Ok(mut service_results) = self.service_results.lock() {
                service_results.push(ServiceResult::Aux(AuxiliaryParameters {
                    energy_standby,
                    energy_crankcase_heater_mode,
                    energy_off_mode,
                }));
            }
            detailed_results.push(self.service_results.clone());
        }

        // Variables below need to be reset at the end of each timestep.
        self.total_time_running_current_timestep = Default::default();
        self.service_results = Default::default();
    }

    pub fn output_detailed_results(
        &self,
        _hot_water_energy_output: f64,
    ) -> (ResultsPerTimestep, ResultsAnnual) {
        let detailed_results = self.detailed_results.as_ref().expect(
            "Detailed results cannot be output when the option to collect them was not selected",
        );
        let mut results_per_timestep: ResultsPerTimestep = Default::default();
        results_per_timestep.insert("auxiliary", Default::default());
        for (parameter, param_unit, _) in AUX_PARAMETERS {
            let mut auxiliary_param_results = results_per_timestep
                .get_mut("auxiliary")
                .unwrap()
                .entry((parameter, param_unit))
                .or_default();
            for service_results in detailed_results {
                if let Ok(service_results) = service_results.lock() {
                    if let ServiceResult::Full(calc) = service_results
                        .iter()
                        .last()
                        .expect("Expected at least one service result")
                    {
                        let result = calc.param(parameter);
                        auxiliary_param_results.push(result);
                    }
                }
            }
        }

        // For each service, report required output parameters
        // TODO iterate over energy supply connections outputting data into results_per_timestep

        let mut results_annual: ResultsAnnual = Default::default();
        results_annual.insert(
            "Overall",
            OUTPUT_PARAMETERS
                .iter()
                .filter(|param| param.2)
                .map(|(parameter, param_units, _)| {
                    (
                        (*parameter, param_units.unwrap_or("")),
                        ResultParamValue::Number(0.0),
                    )
                })
                .collect::<ResultAnnual>(),
        );
        results_annual.insert("auxiliary", Default::default());
        results_annual.get_mut("Overall").unwrap().insert(
            ("energy_delivered_H4", "kWh"),
            ResultParamValue::Number(0.0),
        );
        // Report auxiliary parameters (not specific to a service)
        {
            // make scope for this mutable borrow
            let auxiliary_annual_results = results_annual.get_mut("auxiliary").unwrap();
            for (parameter, param_unit, incl_in_annual) in AUX_PARAMETERS {
                if incl_in_annual {
                    auxiliary_annual_results.insert(
                        (parameter, param_unit),
                        results_per_timestep["auxiliary"][&(parameter, param_unit)]
                            .iter()
                            .copied()
                            .sum::<ResultParamValue>(),
                    );
                }
            }
        }

        // For each service, report required output parameters
        // TODO iterate over energy supply connections outputting data into results_annual

        // Calculate overall CoP for all services combined
        self.calc_service_cop(&mut results_annual);

        (results_per_timestep, results_annual)
    }

    /// Calculate CoP for whole simulation period for the given service (or overall)
    fn calc_service_cop(&self, results_annual: &mut ResultsAnnual) {
        // Add auxiliary energy to overall CoP
        let energy_auxiliary = {
            let results_auxiliary = results_annual.get("auxiliary");
            match results_auxiliary {
                Some(results_aux) => results_aux.values().copied().sum(),
                None => ResultParamValue::Number(0.),
            }
        };

        let mut results_totals = results_annual.get_mut("Overall").unwrap();

        // Calculate CoP at different system boundaries
        let cop_h1_numerator = results_totals[&("energy_delivered_HP", "kWh")];
        let cop_h1_denominator = results_totals[&("energy_input_HP", "kWh")] + energy_auxiliary;
        let cop_h2_numerator = cop_h1_numerator;
        let cop_h2_denominator =
            cop_h1_denominator + results_totals[&("energy_source_circ_pump", "kWh")];
        let cop_h3_numerator =
            cop_h2_numerator + results_totals[&("energy_delivered_backup", "kWh")];
        let cop_h3_denominator =
            cop_h2_denominator + results_totals[&("energy_input_backup", "kWh")];
        let cop_h4_numerator = results_totals[&("energy_delivered_H4", "kWh")];
        let cop_h4_denominator =
            cop_h3_denominator + results_totals[&("energy_heating_circ_pump", "kWh")];

        let cop_h4_note =
            "Note: For water heating services, only valid when HP is only heat source";

        results_totals.insert(("CoP (H1)", ""), cop_h1_numerator / cop_h1_denominator);
        results_totals.insert(("CoP (H2)", ""), cop_h2_numerator / cop_h2_denominator);
        results_totals.insert(("CoP (H3)", ""), cop_h3_numerator / cop_h3_denominator);
        results_totals.insert(
            ("CoP (H4)", cop_h4_note),
            cop_h4_numerator / cop_h4_denominator,
        );
    }
}

const OUTPUT_PARAMETERS: [(&str, Option<&str>, bool); 19] = [
    ("service_name", None, false),
    ("service_type", None, false),
    ("service_on", None, false),
    ("energy_output_required", Some("kWh"), true),
    ("temp_output", Some("K"), false),
    ("temp_source", Some("K"), false),
    ("thermal_capacity_op_cond", Some("kW"), false),
    ("cop_op_cond", None, false),
    ("time_running", Some("hours"), true),
    ("load_ratio", None, false),
    ("hp_operating_in_onoff_mode", None, false),
    ("energy_delivered_HP", Some("kWh"), true),
    ("energy_delivered_backup", Some("kWh"), true),
    ("energy_delivered_total", Some("kWh"), true),
    ("energy_input_HP", Some("kWh"), true),
    ("energy_input_backup", Some("kWh"), true),
    ("energy_heating_circ_pump", Some("kWh"), true),
    ("energy_source_circ_pump", Some("kWh"), true),
    ("energy_input_total", Some("kWh"), true),
];
const AUX_PARAMETERS: [(&str, &str, bool); 3] = [
    ("energy_standby", "kWh", true),
    ("energy_crankcase_heater_mode", "kWh", true),
    ("energy_off_mode", "kWh", true),
];

pub type ResultsPerTimestep<'a> =
    HashMap<&'a str, HashMap<(&'a str, &'a str), Vec<ResultParamValue>>>;
pub type ResultsAnnual<'a> = HashMap<&'a str, ResultAnnual<'a>>;
type ResultAnnual<'a> = HashMap<(&'a str, &'a str), ResultParamValue>;

pub type ResultString = ArrayString<32>;

fn result_str(string: &str) -> ResultString {
    let mut name = ResultString::new();
    name.push_str(string);
    name
}

enum TempSpreadCorrectionArg {
    Float(f64),
    Callable(Box<dyn FnOnce(f64, f64) -> f64>),
}

enum ServiceResult {
    Full(HeatPumpEnergyCalculation),
    Aux(AuxiliaryParameters),
}

struct HeatPumpEnergyCalculation {
    service_name: ResultString,
    service_type: ServiceType,
    service_on: bool,
    energy_output_required: f64,
    temp_output: f64,
    temp_source: f64,
    cop_op_cond: f64,
    thermal_capacity_op_cond: f64,
    time_running: f64,
    deg_coeff_op_cond: f64,
    compressor_power_min_load: f64,
    load_ratio_continuous_min: f64,
    load_ratio: f64,
    use_backup_heater_only: bool,
    hp_operating_in_onoff_mode: bool,
    energy_input_hp_divisor: Option<f64>,
    energy_input_hp: f64,
    energy_delivered_hp: f64,
    energy_input_backup: f64,
    energy_delivered_backup: f64,
    energy_input_total: f64,
    energy_delivered_total: f64,
    energy_heating_circ_pump: f64,
    energy_source_circ_pump: f64,
}

impl HeatPumpEnergyCalculation {
    pub fn param(&self, param: &str) -> ResultParamValue {
        match param {
            "service_name" => ResultParamValue::String(self.service_name),
            "service_type" => ResultParamValue::String(result_str(
                serde_json::to_string(&self.service_type)
                    .unwrap()
                    .as_str()
                    .trim_matches('"'),
            )),
            "service_on" => ResultParamValue::Boolean(self.service_on),
            "energy_output_required" => ResultParamValue::Number(self.energy_output_required),
            "temp_output" => ResultParamValue::Number(self.temp_output),
            "temp_source" => ResultParamValue::Number(self.temp_source),
            "thermal_capacity_op_cond" => ResultParamValue::Number(self.thermal_capacity_op_cond),
            "cop_op_cond" => ResultParamValue::Number(self.cop_op_cond),
            "time_running" => ResultParamValue::Number(self.time_running),
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
            "energy_input_total" => ResultParamValue::Number(self.energy_input_total),
            &_ => panic!("Parameter {param} not recognised"),
        }
    }
}

#[derive(Copy, Clone)]
pub enum ResultParamValue {
    String(ResultString),
    Number(f64),
    Boolean(bool),
}

impl ResultParamValue {
    pub fn as_f64(&self) -> f64 {
        match self {
            ResultParamValue::Number(num) => *num,
            _ => 0.,
        }
    }
}

impl Sum for ResultParamValue {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        // expect everything to be a number (or infer zero), return a number variant
        Self::Number(
            iter.map(|value| match value {
                ResultParamValue::Number(num) => num,
                _ => 0.,
            })
            .sum::<f64>(),
        )
    }
}

impl Add for ResultParamValue {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Number(self.as_f64() + rhs.as_f64())
    }
}

impl Div for ResultParamValue {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self::Number(self.as_f64() / rhs.as_f64())
    }
}

struct AuxiliaryParameters {
    energy_standby: f64,
    energy_crankcase_heater_mode: f64,
    energy_off_mode: f64,
}

/// An object to represent an electric hot-water-only heat pump, tested to EN 16147
#[derive(Clone)]
pub struct HeatPumpHotWaterOnly {
    power_in_kw: f64,
    simulation_timestep: f64,
    control: Option<Arc<Control>>,
    efficiency: f64,
}

impl HeatPumpHotWaterOnly {
    pub fn new(
        power_max_in_kw: f64,
        test_data: &HeatPumpHotWaterTestData,
        vol_daily_average: f64,
        simulation_timestep: f64,
        control: Option<Arc<Control>>,
    ) -> Self {
        let efficiencies = Efficiencies {
            l: test_data.l.as_ref().map(|profile_data| {
                let HeatPumpHotWaterOnlyTestDatum {
                    cop_dhw,
                    hw_tapping_prof_daily,
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
                    hw_tapping_prof_daily,
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
            simulation_timestep,
            control,
            efficiency: Self::init_efficiency(&efficiencies, vol_daily_average),
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

    /// Demand energy (in kWh) from the heat pump
    pub fn demand_energy(&mut self, energy_demand: f64, timestep_idx: usize) -> f64 {
        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        let energy_supplied =
            if self.control.is_none() || self.control.as_ref().unwrap().is_on(timestep_idx) {
                min_of_2(energy_demand, self.power_in_kw * self.simulation_timestep)
            } else {
                0.0
            };

        let _energy_required = energy_supplied / self.efficiency;
        // TODO report to energy supply

        energy_demand
    }

    /// Calculate the maximum energy output (in kWh) from the heater
    pub fn energy_output_max(&self, timestep_idx: usize) -> f64 {
        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        if self.control.is_none() || self.control.as_ref().unwrap().is_on(timestep_idx) {
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
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    pub fn should_interpolate_exhaust_air_heat_pump_test_data() {
        let data_eahp = vec![
            HeatPumpTestDatum {
                air_flow_rate: Some(100.0),
                test_letter: test_letter("A"),
                capacity: 5.0,
                cop: 2.0,
                degradation_coefficient: 0.9,
                design_flow_temp: 55.,
                temp_outlet: 55.,
                temp_source: 20.,
                temp_test: -7.,
            },
            HeatPumpTestDatum {
                air_flow_rate: Some(200.0),
                test_letter: test_letter("A"),
                capacity: 6.0,
                cop: 2.5,
                degradation_coefficient: 0.95,
                design_flow_temp: 55.,
                temp_outlet: 55.,
                temp_source: 20.,
                temp_test: -7.,
            },
            HeatPumpTestDatum {
                air_flow_rate: Some(100.0),
                test_letter: test_letter("B"),
                capacity: 5.5,
                cop: 2.4,
                degradation_coefficient: 0.92,
                design_flow_temp: 35.,
                temp_outlet: 34.,
                temp_source: 20.,
                temp_test: 2.,
            },
            HeatPumpTestDatum {
                air_flow_rate: Some(200.0),
                test_letter: test_letter("B"),
                capacity: 6.0,
                cop: 3.0,
                degradation_coefficient: 0.98,
                design_flow_temp: 35.,
                temp_outlet: 34.,
                temp_source: 20.,
                temp_test: 2.,
            },
        ];
        let data_eahp_interpolated = vec![
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("A"),
                capacity: 5.4,
                cop: 2.2,
                degradation_coefficient: 0.92,
                design_flow_temp: 55.,
                temp_outlet: 55.,
                temp_source: 20.,
                temp_test: -7.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("B"),
                capacity: 5.7,
                cop: 2.64,
                degradation_coefficient: 0.9440000000000001,
                design_flow_temp: 35.,
                temp_outlet: 34.,
                temp_source: 20.,
                temp_test: 2.,
            },
        ];

        let (lowest_air_flow_in_test_data, data_eahp_func_result) =
            interpolate_exhaust_air_heat_pump_test_data(140.0, &data_eahp).unwrap();

        assert_eq!(
            lowest_air_flow_in_test_data, 100.0,
            "incorrect lowest air flow rate identified"
        );
        assert_eq!(
            data_eahp_func_result, data_eahp_interpolated,
            "incorrect interpolation of exhaust air heat pump test data"
        );
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
                test_letter: test_letter("A"),
                capacity: 8.4,
                cop: 4.6,
                degradation_coefficient: 0.90,
                design_flow_temp: 35.,
                temp_outlet: 34.,
                temp_source: 0.,
                temp_test: -7.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("B"),
                capacity: 8.3,
                cop: 4.9,
                degradation_coefficient: 0.90,
                design_flow_temp: 35.,
                temp_outlet: 30.,
                temp_source: 0.,
                temp_test: 2.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("C"),
                capacity: 8.3,
                cop: 5.1,
                degradation_coefficient: 0.90,
                design_flow_temp: 35.,
                temp_outlet: 27.,
                temp_source: 0.,
                temp_test: 7.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("D"),
                capacity: 8.2,
                cop: 5.4,
                degradation_coefficient: 0.95,
                design_flow_temp: 35.,
                temp_outlet: 24.,
                temp_source: 0.,
                temp_test: 12.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("F"),
                capacity: 8.4,
                cop: 4.6,
                degradation_coefficient: 0.90,
                design_flow_temp: 35.,
                temp_outlet: 34.,
                temp_source: 0.,
                temp_test: -7.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("A"),
                capacity: 8.8,
                cop: 3.2,
                degradation_coefficient: 0.90,
                design_flow_temp: 55.,
                temp_outlet: 52.,
                temp_source: 0.,
                temp_test: -7.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("B"),
                capacity: 8.6,
                cop: 3.6,
                degradation_coefficient: 0.90,
                design_flow_temp: 55.,
                temp_outlet: 42.,
                temp_source: 0.,
                temp_test: 2.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("C"),
                capacity: 8.5,
                cop: 3.9,
                degradation_coefficient: 0.98,
                design_flow_temp: 55.,
                temp_outlet: 36.,
                temp_source: 0.,
                temp_test: 7.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("D"),
                capacity: 8.5,
                cop: 4.3,
                degradation_coefficient: 0.98,
                design_flow_temp: 55.,
                temp_outlet: 30.,
                temp_source: 0.,
                temp_test: 12.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("F"),
                capacity: 8.8,
                cop: 3.2,
                degradation_coefficient: 0.90,
                design_flow_temp: 55.,
                temp_outlet: 52.,
                temp_source: 0.,
                temp_test: -7.,
            },
            HeatPumpTestDatum {
                air_flow_rate: None,
                test_letter: test_letter("F2"),
                capacity: 8.8,
                cop: 3.2,
                degradation_coefficient: 0.90,
                design_flow_temp: 55.,
                temp_outlet: 52.,
                temp_source: 0.,
                temp_test: -7.,
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
                    test_letter: test_letter("A"),
                    capacity: 8.4,
                    carnot_cop: 9.033823529411764,
                    cop: 4.6,
                    degradation_coefficient: 0.90,
                    design_flow_temp: 35.,
                    exergetic_eff: 0.5091974605241738,
                    temp_outlet: 34.,
                    temp_source: 0.,
                    temp_test: -7.,
                    theoretical_load_ratio: 1.0,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: test_letter("F"),
                    capacity: 8.4,
                    carnot_cop: 9.033823529438331,
                    cop: 4.6,
                    degradation_coefficient: 0.90,
                    design_flow_temp: 35.,
                    exergetic_eff: 0.5091974605226763,
                    temp_outlet: 34.,
                    temp_source: 0.0000000001,
                    temp_test: -6.9999999999,
                    theoretical_load_ratio: 1.0000000000040385,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: test_letter("B"),
                    capacity: 8.3,
                    carnot_cop: 10.104999999999999,
                    cop: 4.9,
                    degradation_coefficient: 0.90,
                    design_flow_temp: 35.,
                    exergetic_eff: 0.48490846115784275,
                    temp_outlet: 30.,
                    temp_source: 0.,
                    temp_test: 2.,
                    theoretical_load_ratio: 1.1634388356892613,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: test_letter("C"),
                    capacity: 8.3,
                    carnot_cop: 11.116666666666665,
                    cop: 5.1,
                    degradation_coefficient: 0.90,
                    design_flow_temp: 35.,
                    exergetic_eff: 0.4587706146926537,
                    temp_outlet: 27.,
                    temp_source: 0.,
                    temp_test: 7.,
                    theoretical_load_ratio: 1.3186802349509577,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: test_letter("D"),
                    capacity: 8.2,
                    carnot_cop: 12.38125,
                    cop: 5.4,
                    degradation_coefficient: 0.95,
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
                    test_letter: test_letter("A"),
                    capacity: 8.8,
                    carnot_cop: 6.252884615384615,
                    cop: 3.2,
                    degradation_coefficient: 0.90,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.5117638013224666,
                    temp_outlet: 52.,
                    temp_source: 0.,
                    temp_test: -7.,
                    theoretical_load_ratio: 1.0,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: test_letter("F"),
                    capacity: 8.8,
                    carnot_cop: 6.252884615396638,
                    cop: 3.2,
                    degradation_coefficient: 0.90,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.5117638013214826,
                    temp_outlet: 52.,
                    temp_source: 0.0000000001,
                    temp_test: -6.9999999999,
                    theoretical_load_ratio: 1.0000000000030207,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: test_letter("F2"),
                    capacity: 8.8,
                    carnot_cop: 6.252884615408662,
                    cop: 3.2,
                    degradation_coefficient: 0.90,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.5117638013204985,
                    temp_outlet: 52.,
                    temp_source: 0.0000000002,
                    temp_test: -6.9999999998,
                    theoretical_load_ratio: 1.0000000000060418,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: test_letter("B"),
                    capacity: 8.6,
                    carnot_cop: 7.503571428571428,
                    cop: 3.6,
                    degradation_coefficient: 0.90,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.4797715373631604,
                    temp_outlet: 42.,
                    temp_source: 0.,
                    temp_test: 2.,
                    theoretical_load_ratio: 1.3179136223360988,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: test_letter("C"),
                    capacity: 8.5,
                    carnot_cop: 8.587499999999999,
                    cop: 3.9,
                    degradation_coefficient: 0.98,
                    design_flow_temp: 55.,
                    exergetic_eff: 0.4541484716157206,
                    temp_outlet: 36.,
                    temp_source: 0.,
                    temp_test: 7.,
                    theoretical_load_ratio: 1.5978273764295179,
                },
                CompleteHeatPumpTestDatum {
                    air_flow_rate: None,
                    test_letter: test_letter("D"),
                    capacity: 8.5,
                    carnot_cop: 10.104999999999999,
                    cop: 4.3,
                    degradation_coefficient: 0.98,
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
    pub fn should_have_constructed_internal_data_structures(
        test_data: HeatPumpTestData,
        data_sorted: HashMap<OrderedFloat<f64>, Vec<CompleteHeatPumpTestDatum>>,
    ) {
        assert_eq!(
            test_data.dsgn_flow_temps,
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

    #[fixture]
    pub fn expected_init_regression_coeffs() -> HashMap<OrderedFloat<f64>, Vec<f64>> {
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
    pub fn should_populate_init_regression_coeffs(
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
    pub fn should_calc_degradation_coeff(test_data: HeatPumpTestData) {
        let results = [0.9125, 0.919375, 0.92625, 0.933125, 0.94];
        for (i, flow_temp) in [35., 40., 45., 50., 55.].iter().enumerate() {
            assert_eq!(
                round_by_precision(
                    test_data.average_degradation_coeff(celsius_to_kelvin(*flow_temp)),
                    1e7,
                ),
                round_by_precision(results[i], 1e7),
                "incorrect average degradation coefficient returned"
            );
        }
    }

    #[rstest]
    pub fn should_calc_average_capacity(test_data: HeatPumpTestData) {
        let results = [8.3, 8.375, 8.45, 8.525, 8.6];
        for (i, flow_temp) in [35., 40., 45., 50., 55.].iter().enumerate() {
            assert_eq!(
                round_by_precision(
                    test_data.average_capacity(celsius_to_kelvin(*flow_temp)),
                    1e7,
                ),
                round_by_precision(results[i], 1e7),
                "incorrect average capacity returned on step {i}"
            )
        }
    }

    #[rstest]
    pub fn should_calc_temp_spread_test_conditions(test_data: HeatPumpTestData) {
        let results = [5.0, 5.75, 6.5, 7.25, 8.0];
        for (i, flow_temp) in [35., 40., 45., 50., 55.].iter().enumerate() {
            assert_eq!(
                test_data.temp_spread_test_conditions(celsius_to_kelvin(*flow_temp)),
                results[i],
                "incorrect temp spread at test conditions returned"
            );
        }
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
    }

    #[rstest]
    pub fn should_calc_carnot_cop_at_test_condition(
        test_data: HeatPumpTestData,
        carnot_cop_cases: Vec<(f64, &str, f64)>,
    ) {
        for (flow_temp, test_condition, result) in carnot_cop_cases {
            assert_eq!(
                round_by_precision(
                    test_data.carnot_cop_at_test_condition(
                        test_condition,
                        celsius_to_kelvin(flow_temp),
                    ),
                    1e7,
                ),
                round_by_precision(result, 1e7),
                "incorrect Carnot CoP at condition {test_condition} returned"
            );
        }
    }

    #[fixture]
    pub fn outlet_temp_cases() -> Vec<(f64, &'static str, f64)> {
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
    pub fn should_calc_outlet_temp_at_test_condition(
        test_data: HeatPumpTestData,
        outlet_temp_cases: Vec<(f64, &'static str, f64)>,
    ) {
        for (flow_temp, test_condition, result) in outlet_temp_cases {
            assert_eq!(
                round_by_precision(
                    test_data.outlet_temp_at_test_condition(
                        test_condition,
                        celsius_to_kelvin(flow_temp),
                    ),
                    1e7,
                ),
                round_by_precision(result, 1e7),
                "incorrect outlet temp at condition {test_condition} returned"
            );
        }
    }

    #[fixture]
    pub fn source_temp_cases() -> Vec<(f64, &'static str, f64)> {
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
    }

    #[rstest]
    pub fn should_calc_source_temp_at_test_condition(
        test_data: HeatPumpTestData,
        source_temp_cases: Vec<(f64, &'static str, f64)>,
    ) {
        for (flow_temp, test_condition, result) in source_temp_cases {
            assert_eq!(
                round_by_precision(
                    test_data.source_temp_at_test_condition(
                        test_condition,
                        celsius_to_kelvin(flow_temp),
                    ),
                    1e7,
                ),
                round_by_precision(result, 1e7),
                "incorrect source temp at condition {test_condition} returned"
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
    pub fn should_calc_capacity_at_test_condition(
        test_data: HeatPumpTestData,
        capacity_cases: Vec<(f64, &'static str, f64)>,
    ) {
        for (flow_temp, test_condition, result) in capacity_cases {
            assert_eq!(
                round_by_precision(
                    test_data
                        .capacity_at_test_condition(test_condition, celsius_to_kelvin(flow_temp)),
                    1e7,
                ),
                round_by_precision(result, 1e7),
                "incorrect capacity at condition {test_condition} returned"
            );
        }
    }

    #[fixture]
    pub fn lr_op_cond_cases() -> Vec<[f64; 4]> {
        vec![
            [35.0, 283.15, 12.326, 1.50508728516368],
            [40.0, 293.15, 15.6575, 2.38250354792371],
            [45.0, 278.15, 7.95375, 1.21688682087694],
            [50.0, 288.15, 9.23285714285714, 1.58193632324929],
            [55.0, 273.15, 5.96636363636364, 1.0],
        ]
    }

    #[rstest]
    pub fn should_calc_load_ratio_under_operational_conditions(
        test_data: HeatPumpTestData,
        lr_op_cond_cases: Vec<[f64; 4]>,
    ) {
        for [flow_temp, temp_source, carnot_cop_op_cond, result] in lr_op_cond_cases {
            assert_eq!(
                round_by_precision(
                    test_data.load_ratio_at_operating_conditions(
                        celsius_to_kelvin(flow_temp),
                        temp_source,
                        carnot_cop_op_cond,
                    ),
                    1e7,
                ),
                round_by_precision(result, 1e7),
                "incorrect load ratio at operating conditions returned"
            );
        }
    }

    #[rstest]
    pub fn should_calc_load_ratio_degcoeff_either_side_of_op_cond(test_data: HeatPumpTestData) {
        let results_lr_below = [
            1.1634388356892613,
            1.1225791267684564,
            1.0817194178476517,
            1.0408597089268468,
            1.0000000000060418,
            1.3186802349509577,
            1.318488581797243,
            1.3182969286435282,
            1.3181052754898135,
            1.3179136223360988,
        ];
        let results_lr_above = [
            1.3186802349509577,
            1.318488581797243,
            1.3182969286435282,
            1.3181052754898135,
            1.3179136223360988,
            1.513621351820552,
            1.5346728579727933,
            1.555724364125035,
            1.5767758702772765,
            1.5978273764295179,
        ];
        let results_eff_below = [
            0.48490846115784275,
            0.49162229619850667,
            0.49833613123917064,
            0.5050499662798346,
            0.5117638013204985,
            0.4587706146926537,
            0.4640208453602804,
            0.4692710760279071,
            0.4745213066955337,
            0.4797715373631604,
        ];
        let results_eff_above = [
            0.4587706146926537,
            0.4640208453602804,
            0.4692710760279071,
            0.4745213066955337,
            0.4797715373631604,
            0.43614336193841496,
            0.44064463935774134,
            0.4451459167770678,
            0.4496471941963942,
            0.4541484716157206,
        ];
        let results_deg_below = [0.9; 10];
        let results_deg_above = [
            0.9,
            0.9,
            0.9,
            0.9,
            0.9,
            0.95,
            0.9575,
            0.965,
            0.9724999999999999,
            0.98,
        ];

        let mut i = 0;
        for exergy_lr_op_cond in [1.2, 1.4] {
            for flow_temp in [35., 40., 45., 50., 55.] {
                let flow_temp = celsius_to_kelvin(flow_temp);
                let (lr_below, lr_above, eff_below, eff_above, deg_below, deg_above) =
                    test_data.lr_eff_degcoeff_either_side_of_op_cond(flow_temp, exergy_lr_op_cond);
                assert_eq!(
                    round_by_precision(lr_below, 1e7),
                    round_by_precision(results_lr_below[i], 1e7),
                    "incorrect load ratio below operating conditions returned"
                );
                assert_eq!(
                    round_by_precision(lr_above, 1e7),
                    round_by_precision(results_lr_above[i], 1e7),
                    "incorrect load ratio above operating conditions returned"
                );
                assert_eq!(
                    round_by_precision(eff_below, 1e7),
                    round_by_precision(results_eff_below[i], 1e7),
                    "incorrect efficiency below operating conditions returned"
                );
                assert_eq!(
                    round_by_precision(eff_above, 1e7),
                    round_by_precision(results_eff_above[i], 1e7),
                    "incorrect efficiency above operating conditions returned"
                );
                assert_eq!(
                    round_by_precision(deg_below, 1e7),
                    round_by_precision(results_deg_below[i], 1e7),
                    "incorrect degradation coeff below operating conditions returned"
                );
                assert_eq!(
                    round_by_precision(deg_above, 1e7),
                    round_by_precision(results_deg_above[i], 1e7),
                    "incorrect degradation coeff above operating conditions returned"
                );
                i += 1;
            }
        }
    }

    #[rstest]
    pub fn should_calc_cop_op_cond_if_not_air_source(test_data: HeatPumpTestData) {
        let results = [
            6.5629213163133,
            8.09149749487405,
            4.60977003063163,
            5.92554693808559,
            3.76414827675397,
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
            assert_eq!(
                round_by_precision(
                    test_data.cop_op_cond_if_not_air_source(
                        *temp_diff_limit_low,
                        celsius_to_kelvin(*temp_ext),
                        *temp_source,
                        *temp_output,
                    ),
                    1e7,
                ),
                round_by_precision(results[i], 1e7),
                "incorrect CoP at operating conditions (not air source) returned"
            );
        }
    }

    #[rstest]
    pub fn should_calc_capacity_op_cond_if_not_air_source(test_data: HeatPumpTestData) {
        let results = [
            9.26595980986965,
            8.18090909090909,
            8.95014809894768,
            10.0098208201822,
            8.84090909090909,
        ];
        let test_cases = [
            (true, 283.15, 308.15),
            (false, 293.15, 313.15),
            (true, 278.15, 318.15),
            (true, 288.15, 323.15),
            (false, 273.15, 328.15),
        ];
        for (i, (mod_ctrl, temp_source, temp_output)) in test_cases.iter().enumerate() {
            assert_eq!(
                round_by_precision(
                    test_data.capacity_op_cond_if_not_air_source(
                        *temp_output,
                        *temp_source,
                        *mod_ctrl,
                    ),
                    1e7,
                ),
                round_by_precision(results[i], 1e7),
                "incorrect capacity at operating conditions (not air source) returned"
            );
        }
    }

    #[rstest]
    pub fn should_calc_temp_spread_correction(test_data: HeatPumpTestData) {
        let results = [
            1.1219512195122,
            1.08394607843137,
            1.05822498586772,
            1.03966445733223,
            1.02564102564103,
        ];
        let temp_source = 275.15;
        let temp_diff_evaporator = -15.0;
        let temp_diff_condenser = 5.0;
        let temp_spread_emitter = 10.0;

        for (i, temp_output) in [308.15, 313.15, 318.15, 323.15, 328.15].iter().enumerate() {
            assert_eq!(
                round_by_precision(
                    test_data.temp_spread_correction(
                        temp_source,
                        *temp_output,
                        temp_diff_evaporator,
                        temp_diff_condenser,
                        temp_spread_emitter,
                    ),
                    1e7,
                ),
                round_by_precision(results[i], 1e7),
                "incorrect temperature spread correction factor returned"
            );
        }
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }

    fn round_each_by_precision(src: &Vec<f64>, precision: f64) -> Vec<f64> {
        src.iter()
            .map(|num| round_by_precision(*num, precision))
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap()
    }

    impl CompleteHeatPumpTestDatum {
        pub fn round_by_precision(&self, precision: f64) -> CompleteHeatPumpTestDatum {
            CompleteHeatPumpTestDatum {
                air_flow_rate: self.air_flow_rate,
                test_letter: self.test_letter,
                capacity: self.capacity,
                cop: self.cop,
                degradation_coefficient: self.degradation_coefficient,
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

    #[rstest]
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
}
