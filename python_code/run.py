#!/usr/bin/env python3
import glob, argparse, os, json
import numpy as np
from localizer_alignment import localizer_alignment_anat_update_osprey
    
    
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
        registration_output_foler = os.path.join(output_folder, 'PreOspreyLocalizerReg')
        anat_files_dict = localizer_alignment_anat_update_osprey(anat_files_dict, registration_output_foler, localizers)
        
    #Update dictionary with settings for osprey
    joint_dict = settings_dict.copy()
    joint_dict['basisSet'] = os.getenv("BASIS_SETS_PATH")
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

def find_closest_localizer_pair(localizer_pairs_dict, mrs_files_combo, required_suid = None):
    '''Utility that picks which localizer to use with MRS images
    
    This function takes a pair of dictionaries representing localizer and MRS images
    and returns the localizer that is closest in series number to the MRS image in
    the case where the mrs image has a SeriesNumber field in its json sidecar. Otherwise,
    it returns the last localizer(s) in the localizer_pairs_dict.

    localizer_pairs_dict: dictionary of localizer images, where the keys are the series numbers
    and the values are the paths to the images
    mrs_files_combo: dictionary of MRS images, where the values are the paths to the images

    required_suid: None (default, or string). If string, the function will only return localizer pairs that have the same StudyInstanceUID as the MRS images

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
            print('Warning: SeriesNumber field not found associated with MRS file {}, assuming that the last localizer series should be used in registration.\n'.format(mrs_files_combo[mrs_file_req]))

    study_instance_uids = []
    for temp_loc in localizer_pairs_dict.keys():
        localizer_series_numbers.append(temp_loc)
        study_instance_uids.append(nifti_path_to_json_dict(localizer_pairs_dict[temp_loc][0])['StudyInstanceUID'])

    mrs_series_numbers.sort()
    if series_found == 0:
        best_localizer = -1
        localizer_found = False
        if type(required_suid) == type(None):
            for temp_localizer_series in localizer_series_numbers:
                if temp_localizer_series > best_localizer:
                    best_localizer = temp_localizer_series
            return localizer_pairs_dict[best_localizer]
        else:
            if required_suid == "None":
                print('Warning: MRS StudyInstanceUID is set to string "None". Because of this, processing will assume the MRS acquisition and localizer in a given BIDS session directory come from the same scanning session.\n')
            for i, temp_localizer_series in enumerate(localizer_series_numbers):
                if (temp_localizer_series > best_localizer) and ((required_suid == study_instance_uids[i]) or (required_suid == "None")):
                    localizer_found = True
                    best_localizer = temp_localizer_series
            if localizer_found == False:
                raise NameError('Error: no localzier found with same StudyInstanceUID as MRS images. Please check your data. Localizers: {}'.format(localizer_pairs_dict))
            return localizer_pairs_dict[best_localizer]
        
    else:
        best_localizer = None
        localizer_found = False
        print('Warning - script is currently assuming all localizers from a given BIDS imaging session come from the same scanning session.')
        if type(required_suid) == type(None):
            for temp_localizer_series in localizer_series_numbers:
                if temp_localizer_series < mrs_series_numbers[0]:
                    localizer_found = True
                    best_localizer = localizer_pairs_dict[temp_localizer_series]
            if localizer_found == False:
                raise ValueError('Error: Could not find a localizer series that was acquired before the MRS series. Please check your data.')
            else:
                return best_localizer
        else:
            if required_suid == "None":
                print('Warning: MRS StudyInstanceUID is set to string "None". Because of this, processing will assume the MRS acquisition and localizer in a given BIDS session directory come from the same scanning session.\n')
            for i, temp_localizer_series in enumerate(localizer_series_numbers):
                if (temp_localizer_series < mrs_series_numbers[0]) and ((required_suid == study_instance_uids[i]) or (required_suid == "None")):
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


def build_parser():

    #Configure the commands that can be fed to the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("bids_dir", help="The path to the BIDS directory for your study (this is the same for all subjects)", type=str)
    parser.add_argument("output_dir", help="The path to the folder where outputs will be stored (this is the same for all subjects)", type=str)
    parser.add_argument("analysis_level", help="Should always be participant", type=str)
    parser.add_argument("json_settings", help="The path to the subject-agnostic JSON file that will be used to configure processing settings", type=str)

    parser.add_argument('--participant_label', '--participant-label', help="The name/label of the subject to be processed (i.e. sub-01 or 01)", type=str)
    parser.add_argument('--session_id', '--session-id', help="OPTIONAL: the name of a specific session to be processed (i.e. ses-01)", type=str)
    parser.add_argument('--segmentation_dir', '--segmentation-dir', help="The path to the folder where segmentations are stored (this is the same for all subjects)", type=str)
    parser.add_argument('--localizer_registration', '--localizer-registration', help="OPTIONAL: Use localizer to register anatomical images to MRS scan. Also requires the use of --segmentation_dir argument", action='store_true')
    parser.add_argument('--localizer_search_term', '--localizer-search-term', help="OPTIONAL: The search term to use to find localizer images (i.e. *localizer*)", type=str, default='*localizer*.nii*')
    parser.add_argument('--preferred_anat_modality', '--preferred-anat-modality', help="OPTIONAL: The preferred modality to use for anatomical images (i.e. T1w or T2w)", type=str, default='T2w')
    parser.add_argument('--terms_not_allowed_in_anat', nargs='+', help='One or more terms (seperated by spaces) that are not allowed in the file name of high-res anatomical reference images. Useful for getting rid of localizer images that can be confused for high-res anatomical scans. example - "--terms_not_allowed_in_anat mrs ax coronal"', required=False)
    parser.add_argument('--require_same_mrs_localizer_suid', help='If activated, OSPREY will only try to match MRS scans with localizers if the StudyInstanceUID field in the BIDS JSON matches across the localizer/mrs images. Note: if the StudyInstanceUID is undefined for the MRS acquisition, denoted by a value of None, then a localizer will be chosen without taking the StudyInstanceUID into account.', action='store_true')
    
    return parser

def main():

    parser = build_parser()
    args = parser.parse_args()

    compiled_executable_path = os.getenv("EXECUTABLE_PATH")
    mcr_path = os.getenv("MCR_PATH")


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

    if args.preferred_anat_modality not in ['T1w', 'T2w']:
        raise ValueError('Error: preferred_anat_modality must be T1w or T2w, but program received: ' + args.preferred_anat_modality)

    if args.require_same_mrs_localizer_suid:
        require_suid = True
    else:
        require_suid = False

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

            #Grab T1w/T2w file, but only consider files that dont have terms included
            #in args.terms_not_allowed_in_anat
            anats_dict = {}
            t1w_anats_tentative = glob.glob(os.path.join(session_path,'anat/*T1w.ni*'))
            t1w_anats = []
            for temp_t1w in t1w_anats_tentative:
                if args.terms_not_allowed_in_anat:
                    if any(x in temp_t1w.split('/')[-1] for x in args.terms_not_allowed_in_anat):
                        pass
                    else:
                        t1w_anats.append(temp_t1w)
                else:
                    t1w_anats.append(temp_t1w)

            t2w_anats_tentative = glob.glob(os.path.join(session_path,'anat/*T2w.ni*'))
            t2w_anats = []
            for temp_t2w in t2w_anats_tentative:
                if args.terms_not_allowed_in_anat:
                    if any(x in temp_t2w.split('/')[-1] for x in args.terms_not_allowed_in_anat):
                        pass
                    else:
                        t2w_anats.append(temp_t2w)
                else:
                    t2w_anats.append(temp_t2w)

            
            print('Preferred anatomical modality is: ' + args.preferred_anat_modality)
            if (len(t1w_anats) + len(t2w_anats)) == 0:
                print('No T1w or T2w image found for ' + session_path + ', skipping processing for current session.\n')
                continue

            #If T1w is the preferred anatomical reference modality
            elif args.preferred_anat_modality == 'T1w':
                if len(t1w_anats) == 1:
                    anats_dict['files_nii'] = t1w_anats
                elif len(t1w_anats) > 1:
                    raise ValueError('Error: more than 1 T1w image found for ' + session_path + ', processing is not designed to handle use case where multiple T1s are present.')
                elif len(t2w_anats) == 1:
                    anats_dict['files_nii'] = t2w_anats
                elif len(t2w_anats) > 1:
                    raise ValueError('Error: more than 1 T2w image found for ' + session_path + ', processing is not designed to handle use case where multiple T2s are present.')
                else:
                    print('No T1w and T2w image found for ' + session_path + ', skipping processing for current session.')
            
            #If T2w is the preferred anatomical reference modality
            else:
                if len(t2w_anats) == 1:
                    anats_dict['files_nii'] = t2w_anats
                elif len(t2w_anats) > 1:
                    raise ValueError('Error: more than 1 T2w image found for ' + session_path + ', processing is not designed to handle use case where multiple T2s are present.')
                elif len(t1w_anats) == 1:
                    anats_dict['files_nii'] = t1w_anats
                elif len(t1w_anats) > 1:
                    raise ValueError('Error: more than 1 T1w image found for ' + session_path + ', processing is not designed to handle use case where multiple T1s are present.')
                else:
                    print('No T1w and T2w image found for ' + session_path + ', skipping processing for current session.')

            #Grab segmentation
            if isinstance(segmentation_dir, type(None)) == False:
                
                #First look for CABINET output in T2 space

                if 'T2w.nii' in anats_dict['files_nii'][0]:
                    chosen_modality = 'T2w'
                elif 'T1w.nii' in anats_dict['files_nii'][0]:
                    chosen_modality = 'T1w'
                else:
                    raise ValueError('Error: expected to find either T1w or T2w in anatomical file name but found neither. ({})'.format(anats_dict['files_nii'][0]))

                if chosen_modality == 'T2w':
                    t2_seg_path_template = os.path.join(segmentation_dir, temp_participant, temp_session, 'anat', '*_space-T2w_desc-aseg_dseg.nii.gz')
                    t2_seg_files = glob.glob(t2_seg_path_template)
                    if len(t2_seg_files) == 1:
                        anats_dict['files_seg'] = [t2_seg_files]

                if chosen_modality == 'T1w':
                    t1_seg_path_template = os.path.join(segmentation_dir, temp_participant, temp_session, 'anat', '*_space-T1w_desc-aseg_dseg.nii.gz')
                    t1_seg_files = glob.glob(t1_seg_path_template)
                    if len(t1_seg_files) == 1:
                        anats_dict['files_seg'] = [t1_seg_files]

                if 'files_seg' not in anats_dict.keys():
                    raise ValueError('Error: segmentation dir ({}) was provided, but no segmentation was found for subject {} and session {} in either T1w or T2w space.'.format(segmentation_dir, temp_participant, temp_session))
                
                possible_json_file = anats_dict['files_seg'][0][0].split('.nii')[0] + '.json'
                if os.path.exists(possible_json_file):
                    with open(possible_json_file, 'r') as f:
                        json_dict = json.loads(f.read())
                    if 'SpatialReference' in json_dict.keys():
                        if anats_dict['files_nii'][0].split('/')[-1] not in json_dict['SpatialReference']:
                            raise ValueError('Error: segmentation file found for subject {} and session {}, but the SpatialReference field in the json sidecar does not point to the anatomical image chosen for processing.\n\n Segmentation spatial reference: {}\n Reference chosen for current processing {}\n Please check your data.'.format(temp_participant, temp_session, json_dict['SpatialReference'].split('/')[-1], anats_dict['files_nii'][0].split('/')[-1]))
                        else:
                            print('Same spatial reference found for segmentation and anatomical image, proceeding with processing.\n\n Segmentation spatial reference: {}\n Reference chosen for current processing {}\n'.format(json_dict['SpatialReference'].split('/')[-1], anats_dict['files_nii'][0].split('/')[-1]))
                    else:
                        print('No SpatialReference field found in segmentation json sidecar, assuming that segmentation is in the same space as the anatomical image.')
                else:
                    print('No json sidecar found for segmentation, assuming that segmentation is in the same space as the anatomical image.')

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
                print('')

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
                print('File combos: {}'.format(file_combos))
                index = 0
                for temp_combo in file_combos:

                    if require_suid:
                        suid_list = []
                        for temp_mrs_key in temp_combo.keys():
                            temp_mrs_json = nifti_path_to_json_dict(temp_combo[temp_mrs_key][0])
                            try:
                                suid_list.append(temp_mrs_json['StudyInstanceUID'])
                            except:
                                raise ValueError('Error: StudyInstanceUID field not found in BIDS JSONs for MRS images in the current combo. If you still want to process turn argument --require_same_mrs_localizer_suid off')
                        if len(set(suid_list)) > 1:
                            raise ValueError('Error: StudyInstanceUID field in BIDS JSONs for MRS images in the current combo do not match, skipping processing for current combo. File combo: {}'.format(temp_combo))
                        else:
                            current_suid = suid_list[0]
                            print('   Current SUID: {}'.format(current_suid))
                    else:
                        current_suid = None

                    if use_localizer:
                        #Find the best set of localizers for the current MRS images (or use the last localizer)
                        best_localizers = find_closest_localizer_pair(localizer_groups_dict, temp_combo, required_suid = current_suid)
                    else:
                        best_localizers = None

                    #Now need to update run_processing to accept the localizer pair...
                    run_processing(temp_sequence_dict, temp_combo, anats_dict.copy(), output_dir,
                        temp_participant, temp_session, temp_sequence, index,
                        compiled_executable_path, mcr_path, localizers = best_localizers)

                    index += 1

                    



if __name__ == "__main__":
    main()