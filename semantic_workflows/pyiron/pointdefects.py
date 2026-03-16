from pyiron_workflow import as_function_node
from ..pointdefects import (
    create_interstitial as _create_interstitial,
    create_substitutional as _create_substitutional,
    create_vacancy as _create_vacancy,
)
from ..pointdefects import (
    calculate_vacancy_formation_energy,
    calculate_substitutional_formation_energy,
    calculate_interstitial_formation_energy,
)

create_interstitial = as_function_node(_create_interstitial)
create_substitutional = as_function_node(_create_substitutional)
create_vacancy = as_function_node(_create_vacancy)
