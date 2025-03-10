from typing import Dict, Type
from modules.file_conversion.converters.obabel_pdb_to_mol2_converter import (
    OBabelPDBtoMOL2Converter,
)
from modules.file_conversion.converters.base_converter import (
    BaseConverter,
)
from modules.file_conversion.converters.editconf_pdb_to_gro import (
    EditconfPDBtoGROConverter,
)
from modules.file_conversion.converters.editconf_gro_to_pdb import (
    EditconfGROtoPDBConverter,
)


class ConverterFactory:
    _registry: Dict[tuple, Type] = {}

    @classmethod
    def register_converter(cls, converter_class: BaseConverter):
        key = (converter_class.input_file_type, converter_class.output_file_type)
        cls._registry[key] = converter_class

    @classmethod
    def get_converter(cls, input_file_type: str, output_file_type: str):
        key = (input_file_type, output_file_type)
        if key not in cls._registry:
            raise ValueError(f"No converter registered for {key}")
        return cls._registry[key]()


ConverterFactory.register_converter(OBabelPDBtoMOL2Converter)
ConverterFactory.register_converter(EditconfPDBtoGROConverter)
ConverterFactory.register_converter(EditconfGROtoPDBConverter)
