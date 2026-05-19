import os


def update_CARLIN_dir(CARLIN_root_folder, template):
    if template == "cCARLIN":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/cCARLIN'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/cCARLIN"
    elif template == "Tigre":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Tigre_CARLIN'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Tigre_CARLIN"
    elif template == "Tigre_2022":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Tigre_CARLIN_2022'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Tigre_CARLIN_2022"
    elif template == "Tigre_2022_v2":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Tigre_CARLIN_2022_v2'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Tigre_CARLIN_2022_v2"
    elif template == "Rosa":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Rosa_CARLIN'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Rosa_CARLIN"
    elif template == "Rosa_v2":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Rosa_CARLIN_v2'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Rosa_CARLIN_v2"
    else:
        raise ValueError(
            "The input template should be among {Rosa, Tigre_2022, Tigre, cCARLIN, Rosa_v2, Tigre_2022_v2}"
        )
    return os.path.abspath(Actual_CARLIN_dir)
