use fraction::Fraction;

pub struct SimulationTime {
    start_time: Fraction,
    end_time: Fraction,
    step: Fraction,
}

impl SimulationTime {
    fn total_steps(&self) -> Fraction {
        ((self.end_time - self.start_time) / self.step).ceil()
    }

    fn iter(&self) -> SimulationTimeIterator {
        SimulationTimeIterator::from(self)
    }
}

struct SimulationTimeIterator<'a> {
    current_index: usize,
    simulation_time: &'a SimulationTime,
}

impl<'a> SimulationTimeIterator<'a> {
    fn from(simulation_time: &'a SimulationTime) -> Self {
        SimulationTimeIterator {
            current_index: 0,
            simulation_time,
        }
    }
}

struct SimulationTimeIteration {
    pub index: usize,
    pub time: Fraction,
    pub timestep: Fraction,
}

impl<'a> Iterator for SimulationTimeIterator<'a> {
    type Item = SimulationTimeIteration;

    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use rstest::*;
    use fraction::Fraction;
    use super::*;

    #[fixture]
    pub fn timestep() -> Fraction {
        Fraction::from(0.5)
    }

    #[fixture]
    pub fn simtime() -> SimulationTime {
        SimulationTime {
            start_time: Fraction::from(742),
            end_time: Fraction::from(746),
            step: timestep(),
        }
    }

    #[rstest]
    fn should_have_correct_total_steps(simtime: SimulationTime) {
        assert_eq!(simtime.total_steps(), Fraction::from(8))
    }

    #[rstest]
    fn should_iterate_correctly(simtime: SimulationTime, timestep: Fraction) {
        for (i, item) in simtime.iter().enumerate() {
            assert_eq!(item.index, i);
            assert_eq!(item.time, Fraction::from(i) * timestep + Fraction::from(742));
            assert_eq!(item.timestep, timestep);
        }
    }
}
