use crate::compare_floats::max_of_2;
use crate::core::heating_systems::common::SpaceHeatingService;
use crate::core::space_heat_demand::zone::Zone;
use crate::external_conditions::ExternalConditions;
use crate::input::{EcoDesignController, EcoDesignControllerClass};
use crate::simulation_time::SimulationTimeIteration;
// use bacon_sci::ivp::{RungeKuttaSolver, RK45};
use std::sync::Arc;

/// This module provides objects to represent radiator and underfloor emitter systems.

/// Convert flow temperature to return temperature using the 6/7th rule.
///
/// Parameters:
///     `flow_temp_celsius` - Flow temperature in degrees Celsius.
///
///     Returns:
///     float: Return temperature in degrees Celsius.
pub fn convert_flow_to_return_temp(flow_temp_celsius: f64) -> f64 {
    (6.0 / 7.0) * flow_temp_celsius
}

pub struct Emitters {
    thermal_mass: f64,
    c: f64,
    n: f64,
    temp_diff_emit_dsgn: f64,
    frac_convective: f64,
    heat_source: Arc<SpaceHeatingService>,
    zone: Arc<Zone>,
    external_conditions: Arc<ExternalConditions>,
    design_flow_temp: f64,
    ecodesign_controller_class: EcoDesignControllerClass,
    min_outdoor_temp: Option<f64>,
    max_outdoor_temp: Option<f64>,
    min_flow_temp: Option<f64>,
    max_flow_temp: Option<f64>,
    simulation_timestep: f64,
    temp_emitter_prev: f64,
}

impl Emitters {
    /// Construct an Emitters object
    ///
    /// Arguments:
    /// * `thermal_mass` - thermal mass of emitters, in kWh / K
    /// * `c` - constant from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
    /// * `n` - exponent from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
    /// * `temp_diff_emit_dsgn` - design temperature difference across the emitters, in deg C or K
    /// * `frac_convective` - convective fraction for heating
    /// * `heat_source` - reference to an object representing the system (e.g.
    ///                       boiler or heat pump) providing heat to the emitters
    /// * `zone` - reference to the Zone object representing the zone in which the
    ///             emitters are located
    /// * `simulation_timestep` - timestep length for simulation time being used in this context
    ///
    /// Other variables:
    /// * `temp_emitter_prev` - temperature of the emitters at the end of the
    ///    previous timestep, in deg C
    pub fn new(
        thermal_mass: f64,
        c: f64,
        n: f64,
        temp_diff_emit_dsgn: f64,
        frac_convective: f64,
        heat_source: Arc<SpaceHeatingService>,
        zone: Arc<Zone>,
        external_conditions: Arc<ExternalConditions>,
        ecodesign_controller: EcoDesignController,
        design_flow_temp: f64,
        simulation_timestep: f64,
    ) -> Self {
        let ecodesign_controller_class = ecodesign_controller.ecodesign_control_class;
        let (min_outdoor_temp, max_outdoor_temp, min_flow_temp, max_flow_temp) = if matches!(
            ecodesign_controller_class,
            EcoDesignControllerClass::ClassII
                | EcoDesignControllerClass::ClassIII
                | EcoDesignControllerClass::ClassVI
                | EcoDesignControllerClass::ClassVII
        ) {
            (
                ecodesign_controller.min_outdoor_temp,
                ecodesign_controller.max_outdoor_temp,
                ecodesign_controller.min_flow_temp,
                Some(design_flow_temp),
            )
        } else {
            (None, None, None, None)
        };
        Self {
            thermal_mass,
            c,
            n,
            temp_diff_emit_dsgn,
            frac_convective,
            heat_source,
            zone,
            external_conditions,
            design_flow_temp,
            ecodesign_controller_class,
            min_outdoor_temp,
            max_outdoor_temp,
            min_flow_temp,
            max_flow_temp,
            simulation_timestep,
            temp_emitter_prev: 20.0,
        }
    }

    pub fn temp_setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        match self.heat_source.as_ref() {
            SpaceHeatingService::HeatPump(heat_pump) => {
                heat_pump.temp_setpnt(simulation_time_iteration)
            }
            SpaceHeatingService::Boiler(boiler) => boiler.temp_setpnt(*simulation_time_iteration),
            SpaceHeatingService::HeatNetwork(heat_network) => {
                heat_network.temperature_setpnt(simulation_time_iteration)
            }
            SpaceHeatingService::HeatBattery(_) => unreachable!(),
        }
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        match self.heat_source.as_ref() {
            SpaceHeatingService::HeatPump(heat_pump) => {
                heat_pump.in_required_period(simulation_time_iteration)
            }
            SpaceHeatingService::Boiler(boiler) => {
                boiler.in_required_period(*simulation_time_iteration)
            }
            SpaceHeatingService::HeatNetwork(heat_network) => {
                heat_network.in_required_period(simulation_time_iteration)
            }
            SpaceHeatingService::HeatBattery(_) => unreachable!(),
        }
    }

    pub fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    pub fn temp_flow_return(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> (f64, f64) {
        let flow_temp = match self.ecodesign_controller_class {
            EcoDesignControllerClass::ClassII
            | EcoDesignControllerClass::ClassIII
            | EcoDesignControllerClass::ClassVI
            | EcoDesignControllerClass::ClassVII => {
                // A heater flow temperature control that varies the flow temperature of
                // water leaving the heat dependant upon prevailing outside temperature
                // and selected weather compensation curve.
                //
                // They feature provision for manual adjustment of the weather
                // compensation curves and therby introduce a technical risk that optimal
                // minimised flow temperatures are not always achieved.

                // use weather temperature at the timestep
                let outside_temp = self.external_conditions.air_temp(simulation_time_iteration);

                let min_flow_temp = self.min_flow_temp.unwrap();
                let max_flow_temp = self.max_flow_temp.unwrap();
                let min_outdoor_temp = self.min_outdoor_temp.unwrap();
                let max_outdoor_temp = self.max_outdoor_temp.unwrap();

                // set outdoor and flow temp limits for weather compensation curve
                if outside_temp < min_outdoor_temp {
                    max_flow_temp
                } else if outside_temp > max_outdoor_temp {
                    min_flow_temp
                } else {
                    // Interpolate
                    // Note: A previous version used numpy interpolate, but this
                    //        seemed to be giving incorrect results, so interpolation
                    //        is implemented manually here.
                    min_flow_temp
                        + (outside_temp - max_outdoor_temp)
                            * ((max_flow_temp - min_flow_temp)
                                / (min_outdoor_temp - max_outdoor_temp))
                }
            }
            _ => self.design_flow_temp,
        };

        let return_temp = if flow_temp >= 70.0 {
            60.0
        } else {
            flow_temp * 6.0 / 7.0
        };

        (flow_temp, return_temp)
    }

    /// Calculate emitter output at given emitter and room temp
    ///
    /// Power output from emitter (eqn from 2020 ASHRAE Handbook p644):
    ///            power_output = c * (T_E - T_rm) ^ n
    ///        where:
    ///            T_E is mean emitter temperature
    ///            T_rm is air temperature in the room/zone
    ///            c and n are characteristic of the emitters (e.g. derived from BS EN 442 tests)
    pub fn power_output_emitter(&self, temp_emitter: f64, temp_rm: f64) -> f64 {
        self.c * max_of_2(0., temp_emitter - temp_rm).powf(self.n)
    }

    /// Calculate emitter temperature that gives required power output at given room temp
    ///
    /// Power output from emitter (eqn from 2020 ASHRAE Handbook p644):
    ///            power_output = c * (T_E - T_rm) ^ n
    ///        where:
    ///            T_E is mean emitter temperature
    ///            T_rm is air temperature in the room/zone
    ///            c and n are characteristic of the emitters (e.g. derived from BS EN 442 tests)
    ///        Rearrange to solve for T_E
    pub fn temp_emitter_req(&self, power_emitter_req: f64, temp_rm: f64) -> f64 {
        (power_emitter_req / self.c).powf(1. / self.n) + temp_rm
    }

    fn func_temp_emitter_change_rate(
        &self,
        power_input: f64,
    ) -> impl FnOnce(f64, &[f64]) -> f64 + '_ {
        let Self {
            c, n, thermal_mass, ..
        } = self;

        move |t, temp_diff: &[f64]| {
            (power_input - c * max_of_2(0., temp_diff[0]).powf(*n)) / thermal_mass
        }
    }

    // /// Calculate emitter temperature after specified time with specified power input
    // pub fn temp_emitter(
    //     &self,
    //     time_start: f64,
    //     time_end: f64,
    //     temp_emitter_start: f64,
    //     temp_rm: f64,
    //     power_input: f64,
    //     temp_emitter_max: Option<f64>,
    // ) -> Result<(f64, f64), &'static str> {
    //     // Calculate emitter temp at start of timestep
    //
    //     let temp_diff_start = temp_emitter_start - temp_rm;
    // }
    fn demand_energy(&self, energy_demand: f64) -> f64 {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rstest::fixture;
    use rstest::rstest;

    use crate::core::energy_supply::energy_supply::EnergySupply;
    use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
    use crate::core::heating_systems::boiler::Boiler;
    use crate::core::heating_systems::boiler::BoilerServiceWaterRegular;
    use crate::core::heating_systems::common::HeatSourceWet;
    use crate::core::space_heat_demand::zone::Zone;
    use crate::corpus::HeatSource;
    use crate::external_conditions::DaylightSavingsConfig;
    use crate::input::EnergySupplyType;
    use crate::input::HeatSourceWetDetails;
    use crate::simulation_time::SimulationTime;
    use crate::simulation_time::SimulationTimeIterator;

    #[fixture]
    pub(crate) fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0., 2., 0.25).iter()
    }

    #[fixture]
    pub(crate) fn external_conditions() -> ExternalConditions {
        todo!();
    }

    #[fixture]
    pub(crate) fn heat_source() -> SpaceHeatingService {
        todo!();
    }

    #[fixture]
    pub(crate) fn zone() -> Zone {
        todo!();
    }

    #[fixture]
    pub(crate) fn emitters(
        heat_source: SpaceHeatingService,
        zone: Zone,
        external_conditions: ExternalConditions,
    ) -> Emitters {
        let thermal_mass = 0.14;
        let c = 0.08;
        let n = 1.2;
        let temp_diff_emit_dsgn = 10.0;
        let frac_convective = 0.4;

        let ecodesign_controller = EcoDesignController {
            ecodesign_control_class: EcoDesignControllerClass::ClassII,
            min_outdoor_temp: Some(-4.),
            max_outdoor_temp: Some(20.),
            min_flow_temp: Some(30.),
        };

        let design_flow_temp = 55.;

        // TODO check this is correct
        let simulation_timestep = 1.;

        Emitters::new(
            thermal_mass,
            c,
            n,
            temp_diff_emit_dsgn,
            frac_convective,
            heat_source.into(),
            zone.into(),
            external_conditions.into(),
            ecodesign_controller,
            design_flow_temp,
            simulation_timestep,
        )
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    fn test_demand_energy(simulation_time: SimulationTimeIterator, emitters: Emitters) {
        let energy_demand_list = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0];
        let mut energy_demand = 0.0;

        for (t_idx, t_it) in simulation_time.enumerate() {
            energy_demand += energy_demand_list[t_idx];
            let energy_provided = emitters.demand_energy(energy_demand);
            energy_demand -= energy_provided;

            assert_relative_eq!(
                energy_provided,
                [
                    0.26481930394248643,
                    0.8287480680413242,
                    1.053315069769369,
                    1.053315069769369,
                    0.9604801440326911,
                    0.9419772896929609,
                    0.915353814620655,
                    0.7639281136418886
                ][t_idx]
            )
        }
    }
}
