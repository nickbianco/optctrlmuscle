function [] = OpenSim_Muscle_Analysis(motion_file,model_sel,output_path,event)
%OpenSim_Muscle_Analysis Executes a muscle analsysis from the command line
%   Detailed explanation goes here

import org.opensim.modeling.*

% get the analysis tool and change variables
[FunctionPath,~]=fileparts(mfilename('fullpath'));
path_generic_file=fullfile(FunctionPath,'settings_Muscle_analysis.xml');
tool=AnalyzeTool(path_generic_file,false);
tool.setLoadModelAndInput(true)
osimModel=Model(model_sel);
tool.setModel(osimModel);
tool.setResultsDir(output_path);
tool.setInitialTime(event(1));
tool.setFinalTime(event(2));
[~, name, ~]=fileparts(motion_file);
tool.setName(name);

% run the analysis
tool.setModelFilename(model_sel);
tool.setCoordinatesFileName(motion_file);

% if ispc
%     out_path_xml=fullfile(output_path,['muscle_analysis_' name '.xml']);
%     tool.print(out_path_xml);
%     tool.run;
% else
    
    [results_dir,~] = fileparts(output_path);
    curr_dir = fileparts(mfilename('fullpath'));
    if ~exist(fullfile(curr_dir,results_dir), 'dir')
        mkdir(results_dir)
    end
    if ~exist(fullfile(curr_dir,output_path), 'dir')
        mkdir(output_path)
    end
    
    out_path_xml=fullfile(output_path, ['muscle_analysis_' name '.xml']);
    tool.print(out_path_xml);
    %tool_deserialized = AnalyzeTool(out_path_xml);
    %tool_deserialized.run();
    %tool.run;
    %
    % TODOmrs hard-coded OpenSim executable path location.
    assert(~isempty(getenv('OPENSIM_HOME')));
    %assert(exist([getenv('OPENSIM_HOME') '/bin/analyze'], 'file') == 2);
    system(['"' getenv('OPENSIM_HOME') '/bin/analyze" -S ' out_path_xml]);
% end

% dos(['analyze1 -S ', out_path_xml]);

end

