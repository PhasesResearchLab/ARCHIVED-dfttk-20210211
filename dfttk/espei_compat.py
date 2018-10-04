def to_element_case(el):
    """Convert an uppercase elemnt to element case, e.g. FE to Fe, V to V."""
    return el[0].upper() + el[1:].lower()


def dfttk_config_to_espei(config):
    """Convert a DFTTK configuration, e.g. [['Fe', 'Ni'], ['Fe']] to ESPEI's configuration [['FE', 'NI'], 'FE']"""
    espei_config = []
    for subl in config:
        if len(subl) == 1:
            # this sublattice is an endmember
            espei_config.append(subl[0].upper())
        else:
            # this sublattice has an interaction
            espei_config.append([comp.upper() for comp in subl])
    return espei_config


def dfttk_occupancies_to_espei(occupancies):
    """Convert DFTTK occupancies, e.g. [[0.5 0.5], [1.0]] to ESPEI's configuration [[0.5, 0.5], 1.0]"""
    espei_occupancies = []
    for subl in occupancies:
        if len(subl) == 1:
            # this sublattice is an endmember
            espei_occupancies.append(1.0)
        else:
            # this sublattice has an interaction
            espei_occupancies.append([occ for occ in subl])
    return espei_occupancies


def espei_config_to_dfttk(config):
    """Convert a ESPEI configuration, e.g. [['FE', 'NI'], 'FE'] to DFTTK's configuration [['Fe', 'Ni'], ['Fe']]"""
    dfttk_config = []
    for subl in config:
        if isinstance(subl, str):
            # this sublattice is an endmember, e.g. 'FE'
            dfttk_config.append([to_element_case(subl)])
        else:
            # this sublattice is an interacting sublattice, e.g. ['FE', 'NI']
            dfttk_config.append([to_element_case(comp) for comp in subl])
    return dfttk_config


def espei_occupancies_to_dfttk(occupancies):
    """Convert a ESPEI configuration, e.g. [[0.5, 0.5], 1.0] to DFTTK's configuration [[0.5, 0.5], [1.0]]"""

    dfttk_occupancies = []
    for subl in occupancies:
        if isinstance(subl, float):
            # this sublattice is an endmember, e.g. 1.0 occupancy
            dfttk_occupancies.append([1.0])
        else:
            # this sublattice is an interacting sublattice, e.g. [0.5, 0.5]
            dfttk_occupancies.append([f for f in subl])
    return dfttk_occupancies
