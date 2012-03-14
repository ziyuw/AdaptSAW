import ConfigParser

def get_cur_path():
    import os, re
    cur_file = os.path.realpath(__file__)
    cur_path = re.sub('ASAW_config.pyc', '', cur_file)
    cur_path = re.sub('ASAW_config.py', '', cur_path)
    
    return cur_path

def get_lib_dir(config_file=None):
	
	if config_file is None:
		return get_cur_path()
	
	config = ConfigParser.RawConfigParser()
	config.read(get_cur_path() + config_file)
	lib_dir = config.get('Section1', 'lib_dir')
	return lib_dir

def set_cur_res_num(config_file, res_num):
    config = ConfigParser.RawConfigParser()
    config.read(get_cur_path() + config_file)
    if not config.has_section('res_num'):
	config.add_section('res_num')

    config.set('res_num', 'int', str(res_num))
    
    with open(config_file, 'wb') as configfile:
	config.write(configfile)

def get_and_inc_cur_res_num(config_file):
    config = ConfigParser.RawConfigParser()
    config.read(get_cur_path() + config_file)

    if not config.has_section('res_num'):
	set_cur_res_num(config_file, 1)
	config = ConfigParser.RawConfigParser()
	config.read(get_cur_path() + config_file)
	
    cur_res_num = config.getint('res_num', 'int')

    config.set('res_num', 'int', str(cur_res_num+1))
    with open(config_file, 'wb') as configfile:
	config.write(configfile)

    return cur_res_num

def get_cur_res_num(config_file):
    config = ConfigParser.RawConfigParser()
    config.read(get_cur_path() + config_file)

    return config.getint('res_num', 'int')
    
if __name__ == '__main__':
    set_cur_res_num('path_config.cfg', 0)
    #print get_and_inc_cur_res_num('path_config.cfg')