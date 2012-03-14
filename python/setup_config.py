import ConfigParser

config = ConfigParser.RawConfigParser()

home_path = '/home/zwang/project/AdaptSAW/'
config.add_section('Section1')
config.set('Section1', 'lib_dir', home_path+'python/')

with open('path_config.cfg', 'wb') as configfile:
    config.write(configfile)
