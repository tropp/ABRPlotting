import commands

def getcomputer():
    computer_name = commands.getoutput('scutil --get ComputerName')
    if computer_name == 'Lytle':
        basedir = '/Volumes/Pegasus/ManisLab_Data3/abr_data'
    elif computer_name == 'Tamalpais':
        #basedir = '/Users/pbmanis/Desktop/data/ABR_DATA'
        basedir = '/Volumes/Backup2B/ABR_DATA'
    else:
        raise ValueError("Need valid computer name to set base path to data")
    return basedir, computer_name