import os
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TopolGenerator:

    def __init__(self, template_path="modules/gromacs/topol_generator/topol_template"):

        self.template_path = template_path

    def create_topol(self, topol_path: str, res_name: str, itp_path: str):

        if not os.path.exists(self.template_path):
            raise FileNotFoundError(f"Template file not found: {self.template_path}")

        with open(self.template_path, "r") as template_file:
            template_content = template_file.read()

        formatted_content = template_content.format(
            itp_path=f'"{itp_path}"', resname=res_name
        )

        with open(topol_path, "w") as output_file:
            output_file.write(formatted_content)

        logging.info(f"Topology file created: {topol_path}")
        return topol_path
