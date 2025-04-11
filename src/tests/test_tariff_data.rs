mod test_tariff_data {
    use crate::core::energy_supply::tariff_data::*;
    use crate::input::EnergySupplyTariff;
    use rstest::*;
    use std::io::{BufReader, Cursor};

    #[rstest]
    fn test_parse_fixture_file() {
        let data = BufReader::new(Cursor::new(include_str!(
            "../../examples/tariff_data/tariff_data_25-06-2024.csv"
        )));
        let tariff_data = TariffData::new(data).unwrap();
        assert_eq!(
            tariff_data
                .price(&EnergySupplyTariff::VariableTimeOfDay, 10644)
                .unwrap(),
            22.92342657
        );
    }
}
