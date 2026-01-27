from pyiron_workflow import as_function_node
from workflows.pointdefects import (
    create_interstitial as _create_interstitial,
    create_substitutional as _create_substitutional,
)

create_interstitial = as_function_node(_create_interstitial)
create_substitutional = as_function_node(_create_substitutional)
