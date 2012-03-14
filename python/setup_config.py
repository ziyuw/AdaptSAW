import ConfigParser
import ASAW_config

config = ConfigParser.RawConfigParser()

config.add_section('Section1')
config.set('Section1', 'lib_dir', ASAW_config.get_cur_path()+'python/')

with open('path_config.cfg', 'wb') as configfile:
    config.write(configfile)
