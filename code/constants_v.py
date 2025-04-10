"""Constants"""

from os.path import dirname

AMINO_ACIDS = list('ARNDCEQGHILKMFPSTWYV')

ROOT_DIR = dirname(__file__)

LIGATION_DATA = {
    'T4_01h_25C': {
        'file_path':'../data/ligation_data/pmid_30335370/FileS01_T4_01h_25C.csv',
        'enzyme':None,
        'ligase':'T4 Ligase',
        'buffer':'1X T4 DNA Ligase Buffer',
        'incubation_time':'1 hour',
        'incubation_temperature':'25C static',
        'site_size':4,
        'PMID':30335370
    },
    'T4_18h_25C': {
        'file_path':'../data/ligation_data/pmid_30335370/FileS03_T4_18h_25C.csv',
        'enzyme':None,
        'ligase':'T4 Ligase',
        'buffer':'1X T4 DNA Ligase Buffer',
        'incubation_time':'18 hours',
        'incubation_temperature':'25C static',
        'site_size':4,
        'PMID':30335370
    },
    'T4_01h_37C': {
        'file_path':'../data/ligation_data/pmid_30335370/FileS02_T4_01h_37C.csv',
        'enzyme':None,
        'ligase':'T4 Ligase',
        'buffer':'1X T4 DNA Ligase Buffer',
        'incubation_time':'1 hour',
        'incubation_temperature':'37C static',
        'site_size':4,
        'PMID':30335370
    },
    'T4_18h_37C': {
        'file_path':'../data/ligation_data/pmid_30335370/FileS04_T4_18h_37C.csv',
        'enzyme':None,
        'ligase':'T4 Ligase',
        'buffer':'1X T4 DNA Ligase Buffer',
        'incubation_time':'18 hours',
        'incubation_temperature':'37C static',
        'site_size':4,
        'PMID':30335370
    },
    'BsaI_cycling': {
        'file_path':'../data/ligation_data/pmid_32877448/s001_bsai_cycling.csv',
        'enzyme':'BsaI',
        'ligase':'T4 Ligase',
        'buffer':'1X T4 DNA Ligase Buffer',
        'incubation_time':'30 cycles',
        'incubation_temperature':'37C 5 min. -> 16C 5 min. (1 cycle)',
        'site_size':4,
        'PMID':32877448
    },
    'BbsI_cycling': {
        'file_path':'../data/ligation_data/pmid_32877448/s004_bbsi_cycling.csv',
        'enzyme':'BbsI',
        'ligase':'T4 Ligase',
        'buffer':'1X CutSmart Buffer with 10 mM DTT and 1 mM ATP',
        'incubation_time':'30 cycles',
        'incubation_temperature':'37C 5 min. -> 16C 5 min. (1 cycle)',
        'site_size':4,
        'PMID':32877448
    },
    'BsmBI_cycling': {
        'file_path':'../data/ligation_data/pmid_32877448/s002_bsmbi_cycling.csv',
        'enzyme':'BsmBI',
        'ligase':'T4 Ligase',
        'buffer':'1X T4 DNA Ligase Buffer',
        'incubation_time':'30 cycles',
        'incubation_temperature':'42C 5 min. -> 16C 5 min. (1 cycle)',
        'site_size':4,
        'PMID':32877448
    },
    'Esp3I_cycling': {
        'file_path':'../data/ligation_data/pmid_32877448/s003_esp3i_cycling.csv',
        'enzyme':'Esp3I',
        'ligase':'T4 Ligase',
        'buffer':'1X T4 DNA Ligase Buffer',
        'incubation_time':'30 cycles',
        'incubation_temperature':'37C 5 min. -> 16C 5 min. (1 cycle)',
        'site_size':4,
        'PMID':32877448
    },
    'SapI_cycling': {
        'file_path':'../data/ligation_data/pmid_32877448/s005_sapi_cycling.csv',
        'enzyme':'SapI',
        'ligase':'T4 Ligase',
        'buffer':'1X T4 DNA Ligase Buffer',
        'incubation_time':'30 cycles',
        'incubation_temperature':'37C 5 min. -> 16C 5 min. (1 cycle)',
        'site_size':3,
        'PMID':32877448
    }
}
