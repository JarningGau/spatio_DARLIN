import os
help_function_dir, this_filename = os.path.split(__file__)
root_dir=os.path.dirname(help_function_dir)
script_dir=os.path.join(root_dir,'bin')
QC_dir=os.path.join(root_dir,'QC')
ref_dir=os.path.join(root_dir,'reference')
