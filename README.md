# NeuroScienceProject
  Preprocessing algorithm for TEP experiment EEG data.
  
  This code flow locates the TMS pulses contaminated data, removes them and replaces them by cubic interpolation.
  
  Full description of the project's goals in 'Research problem.docx'.

# To run locally
  Open cmdline in your folder and run : pip install -r requirements.txt
  
  Put  'get_pulses', 'implement_interpolation' in your folder and add "from implement_interpolation import *" to the head of your file
  
  call tms_pulse_interpolation(raw_or_epoch_object, indices_before_pulse, indices_after_pulse)

# Files description
  'Area Mapping.txt'
      A list of all row indices with their coresponding brain region
      
  'Research problem.docx'
      A one page description of this project's goals
      
  'draw_all.py'
      A code designed to create folders containing plots for each brain region seperatly:
      
      Complete trial, Pulses zoom out, (Pulse zoom in) * No. of pulses
      
      Input: matrix, list of indices to focus on, folder name, list of channels/segments
      
      Output: Folder for each region with plots
      
  'get_pulses'
      A code flow designated to find the location of the TMS pulses.
      
      Input: Raw object #current version creates raw object from hard-coded file path
      
      Output: A log dictionary containing the indices to interpolate.
      
              The indices between each pair are the indices which their values will be removed and replaced by interpolation
    
   'implement_interpolation.py'
   
    The main flow
    
        
        Calls get_pulses, creates the log, creates new interpolated matrix and plot by request
        
        Output:
          Plot of one example segment of the interpolated EEG data
          A new object with an interpolated matrix

    'external.py' 
      An example file for how we created the objects and ran the functions.
          
 # Links & Contacts:
   
    amitomer312@gmail.com
    
    amithorovitz@mail.tau.ac.il
    
    matanshvide@mail.tau.ac.il
    :)
    
    

