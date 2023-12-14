#!/usr/bin/env python3
import glob, argparse, os, json
import numpy as np
from localizer_alignment import localizer_alignment_anat_update_osprey

#To change
#participant label vs sub - done
#run on all participants - done
#take analysis level - done
#return the boutiques descriptor?
#execution.json?

#Configure the commands that can be fed to the command line
parser = argparse.ArgumentParser()
parser.add_argument("bids_dir", help="The path to the BIDS directory for your study (this is the same for all subjects)", type=str)
parser.add_argument("output_dir", help="The path to the folder where outputs will be stored (this is the same for all subjects)", type=str)
parser.add_argument("analysis_level", help="Should always be participant", type=str)
parser.add_argument("json_settings", help="The path to the subject-agnostic JSON file that will be used to configure processing settings", type=str)

parser.add_argument('--participant_label', '--participant-label', help="The name/label of the subject to be processed (i.e. sub-01 or 01)", type=str)
parser.add_argument('--segmentation_dir', '--segmentation-dir', help="The path to the folder where segmentations are stored (this is the same for all subjects)", type=str)
parser.add_argument('--session_id', '--session-id', help="OPTIONAL: the name of a specific session to be processed (i.e. ses-01)", type=str)
parser.add_argument('--localizer_registration', '--localizer-registration', help="OPTIONAL: Use localizer to register anatomical images to MRS scan. Also requires the use of --segmentation_dir argument", action='store_true')
parser.add_argument('--localizer_search_term', '--localizer-search-term', help="OPTIONAL: The search term to use to find localizer images (i.e. *localizer*)", type=str, default='*localizer*.nii*')
args = parser.parse_args()

compiled_executable_path = os.getenv("EXECUTABLE_PATH")
mcr_path = os.getenv("MCR_PATH")
basis_sets_path = os.getenv("BASIS_SETS_PATH")


#export EXECUTABLE_PATH=/home/midb-ig/shared/HBCD/mrs/osprey_dirs_for_containerization/osprey_containerization_code_v1/run_compiled.sh
#export MCR_PATH=/home/midb-ig/shared/repositories/leex6144/v910

#Get cwd in case relative paths are given
cwd = os.getcwd()

#reassign variables to command line input
bids_dir = args.bids_dir
if os.path.isabs(bids_dir) == False:
	bids_dir = os.path.join(cwd, bids_dir)
output_dir = args.output_dir
if os.path.isabs(output_dir) == False:
	output_dir = os.path.join(cwd, output_dir)
analysis_level = args.analysis_level
if analysis_level != 'participant':
    raise ValueError('Error: analysis level must be participant, but program received: ' + analysis_level)
json_settings = args.json_settings
if os.path.isabs(json_settings) == False:
	json_settings = os.path.join(cwd, json_settings)


#Set session label
if args.session_id:
    session_label = args.session_id
    if 'ses-' not in session_label:
        session_label = 'ses-' + session_label
else:
    session_label = None
    
#Set segmentation directory
if args.segmentation_dir:
    segmentation_dir = args.segmentation_dir
    if os.path.isabs(segmentation_dir) == False:
        segmentation_dir = os.path.join(cwd, segmentation_dir)
    use_localizer = args.localizer_registration
else:
    segmentation_dir = None
    use_localizer = False
    
#Find participants to try running
if args.participant_label:
    participant_split = args.participant_label.split(' ')
    participants = []
    for temp_participant in participant_split:
        if 'sub-' not in temp_participant:
            participants.append('sub-' + temp_participant)
        else:
            participants.append(temp_participant)
else:
    os.chdir(bids_dir)
    participants = glob.glob('sub-*')
    
    
    
#############################################################################################################
#############################################################################################################
#############################################################################################################

def run_processing(settings_dict, mrs_files_dict, anat_files_dict, derivs_folder_path,
                   participant_label, session_partial_path, sequence, index,
                   compiled_executable_path, mcr_path, localizers = None):
        
    #If this is the first instance of this session/config combo being used for processing
    #dont number it with *_index, otherwise add the index
    if index == 0:
        output_folder = os.path.join(derivs_folder_path, participant_label, session_partial_path, sequence)
    else:
        output_folder = os.path.join(derivs_folder_path, participant_label, session_partial_path, sequence + '_filecombo_' + str(index))
        
    #Create output folder to store results and json config file
    if os.path.exists(output_folder) == False:
        os.makedirs(output_folder)

    #If localizer registration is being used, register the anat image + segmentation to
    #the localizer, and update the anat_files_dict to reference the newly registered images
    if type(localizers) == type(None):
        pass
    else:
        registration_output_foler = os.path.join(output_folder, 'pre_osprey_localizer_reg')
        anat_files_dict = localizer_alignment_anat_update_osprey(anat_files_dict, registration_output_foler, localizers)
        
    #Update dictionary with settings for osprey
    joint_dict = settings_dict.copy()
    joint_dict['basisSet'] = basis_sets_path
    joint_dict.update(mrs_files_dict)
    joint_dict.update(anat_files_dict)
    del joint_dict['prerequisites']
    joint_dict['outputFolder'] = [output_folder]
    
    #Save settings to output directory
    json_output_path = os.path.join(output_folder, 'wrapper_settings.json')
    with open(json_output_path, 'w') as f:
        joint_json = json.dumps(joint_dict, indent = 4)
        f.write(joint_json)
    
    #Run osprey
    print('Running: ' + compiled_executable_path + ' ' + mcr_path + ' ' + json_output_path)
    output_status = os.system(compiled_executable_path + ' ' + mcr_path + ' ' + json_output_path)
    output_status = os.WEXITSTATUS(output_status)
    if output_status > 0:
        raise ValueError('Error: Matlab command line returned non-zero exit status ({})'.format(output_status))
    
    return

def nifti_path_to_json_dict(nifti_path):
    
    #Function that takes the path to a nifti
    #and returns the contents of its json
    #sidecar as a dictionary
    
    pre_name = nifti_path.split('.nii')[0]
    json_name = pre_name + '.json'
    if os.path.exists(json_name) == False:
        raise AttributeError('Error: Could not find the expected json file ' + json_name)
    else:
        with open(json_name, 'r') as f:
            json_dict = json.loads(f.read())
            
    return json_dict

def find_closest_localizer_pair(localizer_pairs_dict, mrs_files_combo):
    '''Utility that picks which localizer to use with MRS images
    
    This function takes a pair of dictionaries representing localizer and MRS images
    and returns the localizer that is closest in series number to the MRS image in
    the case where the mrs image has a SeriesNumber field in its json sidecar. Otherwise,
    it returns the last localizer(s) in the localizer_pairs_dict.

    localizer_pairs_dict: dictionary of localizer images, where the keys are the series numbers
    and the values are the paths to the images
    mrs_files_combo: dictionary of MRS images, where the values are the paths to the images

    returns: list of paths to localizer images to use for registration
    
    '''

    mrs_series_numbers = []
    localizer_series_numbers = []
    series_found = 0
    for mrs_file_req in mrs_files_combo.keys():
        try:
            mrs_series_numbers.append(nifti_path_to_json_dict(mrs_files_combo[mrs_file_req])['SeriesNumber'])
            series_found = 1
        except:
            print('Warning: SeriesNumber field not found associated with MRS file {}, assuming that the last localizer series should be used in registration.'.format(mrs_files_combo[mrs_file_req]))

    for temp_loc in localizer_pairs_dict.keys():
        localizer_series_numbers.append(temp_loc)


    mrs_series_numbers.sort()
    localizer_series_numbers.sort()
    if series_found == 0:
        return localizer_pairs_dict[localizer_series_numbers[-1]]
        
    else:
        best_localizer = None
        localizer_found = False
        for temp_localizer_series in localizer_series_numbers:
            if temp_localizer_series < mrs_series_numbers[0]:
                localizer_found = True
                best_localizer = localizer_pairs_dict[temp_localizer_series]
        if localizer_found == False:
            raise ValueError('Error: Could not find a localizer series that was acquired before the MRS series. Please check your data.')
        else:
            return best_localizer

    return


def find_acceptable_file_combos(prereq_dict, subj_dir):
    '''Function that finds files paired by intendedfor
    
    The input to this function is a dictionary whose keys
    represent specific individual requirements, and values represent
    nifti files that fill those requirements. The function
    will then iterate through all combinations of requirements
    to see what permutations exist that fullfill all requirements
    at the same time. In the case where there are any fields without
    files that satisfy the requirements, the function will return an
    empty array. In the case where there is exactly one file per requirement,
    the function will assume that those files can be used together. In
    the case where any requirement are associated with multiple files,
    then the function will (1) look at the json sidecar that accompanies
    the nifti file (i.e. swap out .nii or .nii.gz with .json) and then
    (2) look for the IntendedFor field within the json sidecar, for (3)
    all possible solutions where the requirements for each field have
    IntendedFors that include the requirements for the other fields.
    
    This means that if a file that fills one requirement has an intended for
    that points to a file that fills a second requirement, the second file must
    have an intendedfor that also points to the first file for the pair to be
    used in processing. 
    
    IntendedFor fields should start at the session level (i.e. ses/func/run.nii.gz)
    
    Params
    ------
    
    prereq_dict : dict
        A dictionary whose keys represent requirements, and values are lists
        of files that fullfil those requirements
        
    Returns
    -------
    
    list of dicts where the dicts have one file per requirement, (values are lists of
    strings, but the assumed behavior will be that each list has one file).
    '''
    
    #Figure out how many files are present for
    #each requirement
    key_lens = []
    keys = list(prereq_dict.keys())
    for temp_key in keys:
        key_lens.append(len(prereq_dict[temp_key]))
        
    #If any requirement isn't satisfied then
    #exit
    if any(x == 0 for x in key_lens):
        return []
        
    #If there is one file for each requirement,
    #don't look through intendedfor fields
    if all(x == 1 for x in key_lens):
        
        #Change the prereqs to have full paths instead
        #of being relative to the subject dir
        for temp_key in prereq_dict.keys():
            requirements_temp_list = prereq_dict[temp_key]
            for i in range(len(requirements_temp_list)):
                prereq_dict[temp_key][i] = os.path.join(subj_dir, prereq_dict[temp_key][i])
                
        return [prereq_dict]

    trash = np.zeros(key_lens)
    #This is a list of tuples. For each list index,
    #the first index of the underlying tuple will
    #have the values necessary to iterate through
    #all combinations of files in prereq_dict,
    #choosing 1 file per requirement.
    array_enumerator = list(np.ndenumerate(trash))
    
    acceptable_combinations = []
    for i in array_enumerator:
        
        #This is a tuple that will tell
        #you which files from prereq_dict
        #to grab for the current iteration.
        #With these files, we will then check
        #whether the intendedfors of the various
        #requirements work with one another
        x = i[0]

        intendedfors = []
        nifti_names = []
        temp_dict = {}
        for j, k in enumerate(x):
            temp_nifti = prereq_dict[keys[j]][k]
            temp_json = nifti_path_to_json_dict(temp_nifti)
            
            if 'IntendedFor' in temp_json.keys():
                intendedfors.append(temp_json['IntendedFor'])
                temp_dict[keys[j]] = [os.path.join(subj_dir,temp_nifti)] #THIS IS IMPLEMENTED HERE BUT DOESN"T TAKE EFFECT FOR SINGLE FILE USE CASE (this is okay)
                nifti_names.append(temp_nifti) 
            else:
                raise ValueError('Error: if any requirement has more than one file, IntendedFor fields from the json sidecar must be present to proceed. IntendedFor field was not found for json associated with: ' + temp_nifti)

        process = True
        for j, temp_nifti in enumerate(nifti_names):

            intendedfors_without_current_key = intendedfors.copy()
            intendedfors_without_current_key.pop(j)
            
            #This is for the case where there is only one requirement file
            if len(intendedfors_without_current_key) == 0:
                continue
            #This is for the case where there are multiple requirement files
            else:
                for temp_intendedfor in intendedfors_without_current_key:
                    if temp_nifti in temp_intendedfor:
                        pass
                    else:
                        process = False

        if process:
            acceptable_combinations.append(temp_dict)
            
    return acceptable_combinations


#############################################################################################################
#############################################################################################################
#############################################################################################################

#Load json settings file
with open(json_settings, 'r') as f:
    master_settings = json.loads(f.read())


#Iterate through all participants
for temp_participant in participants:
    
    #Check that participant exists at expected path
    subject_path = os.path.join(bids_dir, temp_participant)
    if os.path.exists(subject_path):
        os.chdir(subject_path)
    else:
        raise AttributeError('Error: no directory found at: ' + subject_path)
    
    #Find session/sessions
    if session_label == None:
        sessions = glob.glob('ses*')
        if len(sessions) < 1:
            sessions = ['']
    elif os.path.exists(session_label):
        sessions = [session_label]
    else:
        raise AttributeError('Error: session with name ' + session_label + ' does not exist at ' + subject_path)

    #Iterate through sessions
    for temp_session in sessions:

        #If there is no session structure, this will go to the subject path
        session_path = os.path.join(subject_path, temp_session)

        #Grab T1w file
        anats_dict = {}
        anats = glob.glob(os.path.join(session_path,'anat/*T1w.ni*'))
        if len(anats) == 0:
            print('No T1w image found for ' + session_path + ', skipping processing for current session.')
            continue
        elif len(anats) == 1:
            anats_dict['files_nii'] = anats
        else:
            raise ValueError('Error: more than 1 T1w image found for ' + session_path + ', processing is not designed to handle use case where multiple T1s are present.')

        #Grab segmentation
        if isinstance(segmentation_dir, type(None)) == False:
            
            #First look for CABINET output in T2 space
            t2_seg_path_template = os.path.join(segmentation_dir, temp_participant, temp_session, 'anat', '*_space-T2w_desc-aseg_dseg.nii.gz')
            t2_seg_files = glob.glob(t2_seg_path_template)
            t1_seg_path_template = os.path.join(segmentation_dir, temp_participant, temp_session, 'anat', '*_space-T1w_desc-aseg_dseg.nii.gz')
            t1_seg_files = glob.glob(t1_seg_path_template)
            space_seg_path_template = os.path.join(segmentation_dir, temp_participant, temp_session, 'anat', '*_space-orig_desc-aseg_dseg.nii.gz')
            space_seg_files = glob.glob(space_seg_path_template)
            if len(t2_seg_files) == 1:
                anats_dict['files_seg'] = [t2_seg_files]
            elif len(t1_seg_files) == 1:
                anats_dict['files_seg'] = [t1_seg_files]
            elif len(space_seg_files) == 1:
                anats_dict['files_seg'] = [space_seg_files]
            else:
                raise ValueError('Error: expected to find exactly one segmentation in either T2w, T1w, or orig space but found: {} {} {}, respectively (at least 1 needs to have 1).'.format(len(t2_seg_files), len(t1_seg_files), len(space_seg_files)))

        #Iterate through processing configurations defined in the input json file
        for temp_sequence in master_settings.keys():

            os.chdir(subject_path)
            temp_sequence_dict = master_settings[temp_sequence]
            prereq_keys = temp_sequence_dict['prerequisites'].keys()
            prereq_dict = {}
            num_files = []
            all_files = []
            print("\nFor session " + session_path + " and configuration " + temp_sequence + ":")
            for temp_prereq in prereq_keys:
                prereq_dict[temp_prereq] = glob.glob(os.path.join(temp_session, 'mrs', temp_sequence_dict['prerequisites'][temp_prereq]))
                print("Found " + str(len(prereq_dict[temp_prereq])) + " " + temp_prereq)
                num_files.append(len(prereq_dict[temp_prereq]))
                all_files += prereq_dict[temp_prereq] 

            #Grab all localizer scans within the current session and organize them into groups based on SeriesNumber
            if use_localizer:
                localizer_imgs = glob.glob(os.path.join(session_path, 'anat/{}'.format(args.localizer_search_term)))
                if len(localizer_imgs) == 0:
                    print('No localizer images found for ' + session_path + ', skipping processing for current session.')
                    continue
                else:
                    localizer_series_nums = []
                    for temp_img in localizer_imgs:
                        temp_json = nifti_path_to_json_dict(temp_img)
                        localizer_series_nums.append(temp_json['SeriesNumber'])
                    unique_localizer_series_nums = list(set(localizer_series_nums))
                    localizer_groups_dict = {}
                    for temp_series_num in unique_localizer_series_nums:
                        localizer_groups_dict[temp_series_num] = []
                    for t, temp_series_num_sub in enumerate(localizer_series_nums):
                        localizer_groups_dict[temp_series_num_sub].append(localizer_imgs[t])



            file_combos = find_acceptable_file_combos(prereq_dict, subject_path)
            index = 0
            for temp_combo in file_combos:


                if use_localizer:
                    #Find the best set of localizers for the current MRS images (or use the last localizer)
                    best_localizers = find_closest_localizer_pair(localizer_groups_dict, temp_combo)
                else:
                    best_localizers = None

                #Now need to update run_processing to accept the localizer pair...

                run_processing(temp_sequence_dict, temp_combo, anats_dict.copy(), output_dir,
                       temp_participant, temp_session, temp_sequence, index,
                       compiled_executable_path, mcr_path, localizers = best_localizers)

                index += 1
