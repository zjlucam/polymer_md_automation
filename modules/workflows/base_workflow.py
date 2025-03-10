from abc import ABC, abstractmethod


class BaseWorkflow(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def run(self, **kwargs) -> str:
        pass
