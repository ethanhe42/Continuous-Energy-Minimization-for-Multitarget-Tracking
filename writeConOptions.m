function status = writeConOptions(opt,inifile)
% write ini config file

ini = IniConfig();
ini.AddSections('Parameters');
ini.AddSections('General');
ini.AddSections('Initialization');
ini.AddSections('Detections');
ini.AddSections('Appearance');
ini.AddSections('Miscellaneous');

opttmp=opt;

while ~isempty(fieldnames(opttmp))
    fnames=fieldnames(opttmp);
    fieldname=fnames{1};
    ini=parseField(ini,opt,fieldname);
    opttmp=rmfield(opttmp,fieldname);
end

status = ini.WriteFile(inifile);

end

function ini=parseField(ini,opt,fieldname)
    fvalue=getfield(opt,fieldname);
    if isstruct(fvalue)
        opttmp=fvalue;
        switch (fieldname)
            % Detections
            case 'detScale'
                while ~isempty(fieldnames(opttmp))
                    fnames=fieldnames(opttmp);
                    fieldname=fnames{1};
                    fvalue=getfield(opttmp,fieldname);
                    ini.AddKeys('Detections',fieldname,fvalue);
                    opttmp=rmfield(opttmp,fieldname);
                end                
            % Appearance Parameters
            case 'app'
                while ~isempty(fieldnames(opttmp))
                    fnames=fieldnames(opttmp);
                    fieldname=fnames{1};
                    fvalue=getfield(opttmp,fieldname);
                    ini.AddKeys('Appearance',fieldname,fvalue);
                    opttmp=rmfield(opttmp,fieldname);
                end                             
        end
        % frames, special case
    elseif strcmpi(fieldname,'frames')
        ini.AddKeys('General','ff',fvalue(1));
        ini.AddKeys('General','lf',fvalue(end));
    elseif strcmpi(fieldname,'seqLength')
        % ignore
    else
        sec=getSectionName(fieldname);
        ini.AddKeys(sec,fieldname,fvalue);
    end    
end

function sec=getSectionName(fieldname)
sec='Miscellaneous';
switch(fieldname)
    % Parameters
    case {'wtEdet', ...
            'wtEdyn', ...
            'wtEexc', ...
            'wtEper', ...
            'wtEreg', ...
            'wtEapp', ...
            'lambda' ...
          }
        sec='Parameters';
        
        % General
    case {'track3d',...
            'verbosity', ...	
            'mex', ...		
            'visOptim', ...	
            'cutToTA', ...	
            'remOcc', ...		
            'maxIterCGD', ...	
            'occ', ...		
            'jumpsOrder', ...
            'maxEpochs' ...
         }
        sec='General';
    case {  'startsol', ...
            'EKFDir', ...
            'DPDir' ...            
            }
        sec='Initialization';
    case {
            'detThreshold', ...
            'sigA', ...
            'sigB' ...
            }
        sec='Detections';
   
end
        

end