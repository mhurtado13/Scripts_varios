# modules
import argparse
import os
import shutil
import sys

# const
program_name = sys.argv[0]
program_name_without_extension = os.path.splitext(program_name)[0]


# function
def install_pcuserproj(model_name, gym_name, s_root='../PhysiCell/', b_force=False):
    """
    input:
        check out argparse in __main__ !

    output:
        PhysiCell/user_projects/model_name/ folder.

    description:
        check out argparse in __main__ !
    """

    # check for physicell root folder
    s_root = s_root.replace('\\','/')
    if not s_root.endswith('/'):
        s_root = s_root + '/'
    if not (os.path.isdir(s_root)):
        sys.exit(f"Error @ install_pcuserproj : no PhysiCell root directory found at '{s_root}'.")

    # check for model project folder
    s_path_prj = f'{s_root}user_projects/{model_name}/'
    if ((not b_force) and os.path.exists(s_path_prj)):
        sys.exit(f"Error @ install_pcuserproj : {s_path_prj} already exists!\nUse the command line -f or --force argument to overwrite the existing model.")

    # install project
    print(f'install {s_path_prj} ...')
    if os.path.exists(s_path_prj):
        shutil.rmtree(s_path_prj)

    # copy the entire folder into PhysiCell user_project
    shutil.copytree(model_name, s_path_prj)
    print(f"{model_name} installed!")

    # check for gym project folder
    s_path_prj = f'{s_root}{gym_name}/'
    if ((not b_force) and os.path.exists(s_path_prj)):
        sys.exit(f"Error @ install_pcuserproj : {s_path_prj} already exists!\nUse the command line -f or --force argument to overwrite the existing RL folder.")

    # install project
    print(f'install {s_path_prj} ...')
    if os.path.exists(s_path_prj):
        shutil.rmtree(s_path_prj)

    # copy the entire folder into PhysiCell user_project
    shutil.copytree(gym_name, s_path_prj)

    print("gym env also installed!")


# run
if __name__ == "__main__":
    print('running install_model.py script ...')

    # argv
    parser = argparse.ArgumentParser(
        prog = 'install ',
        description = 'script to copy the PhysiCell model user project into the correct folder structure.',
        epilog = 'afterwards RL can be built, run, and further developed within PhysiCell as usual.',
    )

    # project_name
    parser.add_argument(
        'model_name',
        help='Name of the project to install.'
    )
    # gym_project_name
    parser.add_argument(
        'gym_name',
        help='Name of the RL project to install.'
    )

    # s_root
    parser.add_argument(
        'root',
        nargs = '?',
        default = '../PhysiCell/',
        help = 'path to the PhysiCell root directory.'
    )
    # b_force
    parser.add_argument(
        '-f', '--force',
        #type = bool,
        #nargs = 0,
        action=argparse.BooleanOptionalAction,
        #default = False,
        help = 'Overwrite the existing model with this RL_TME model.'
    )

    # parse arguments
    args = parser.parse_args()
    #print(args)

    # processing
    install_pcuserproj(args.model_name, args.gym_name, args.root, args.force)
