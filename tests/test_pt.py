#!python
#

def get_density_from_pt(ele_list):
    """
    Get density from periodictable package

    Parameters
    ----------
        ele_list : list
            The list of elements, e.g. ['V', 'Ni']
    Returns
    -------
        density_dict : dict
            Dictionary of {element: density}, e.g. {'V': 6.313, 'Ni': 9.03}. 
    Examples
    --------
    >>> get_density_from_pt(['V', 'Ni'])
    {'V': 6.313, 'Ni': 9.03}
    """
    density_dict = {}
    import periodictable as pt
    for ele in ele_list:
        density_dict[ele] = eval("pt." + ele + ".density")
    return density_dict

print(get_density_from_pt(['Nb', 'Ti']))