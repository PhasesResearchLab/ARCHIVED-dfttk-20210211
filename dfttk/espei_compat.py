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


def espei_config_to_dfttk(config):
    """Convert a ESPEI configuration, e.g. [['FE', 'NI'], 'FE'] to DFTTK's configuration [['Fe', 'Ni'], ['Fe']]"""
    def element_case(el):
        """Convert an uppercase elemnt to element case, e.g. FE to Fe, V to V."""
        return el[0].upper() + el[1:].lower()

    dfttk_config = []
    for subl in config:
        if isinstance(subl, str):
            # this sublattice is an endmember, e.g. 'FE'
            dfttk_config.append([element_case(subl)])
        else:
            # this sublattice is an interacting sublattice, e.g. ['FE', 'NI']
            dfttk_config.append([element_case(comp) for comp in subl])
    return dfttk_config
