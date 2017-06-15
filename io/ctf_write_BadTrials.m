function [newexcludelength] = ctf_write_BadTrials(exclude, Dataset, maxexclude)

% This should write a list of badtrials directly into the ClassFile.cls.
% 
% USAGE: [newexcludelength] = ctf_write_BadTrials(exclude, Dataset, maxexclude)
%
% INPUTS:
%   exclude = a cell array (1 x number of trials) listing trials to exclude
%   Dataset = full path of the dataset in question
%   maxexclude = the maximum number of trials to exclude before reverting
%   to only the existing marked trials. Use 50 if unsure, it's just a
%   safety thing.
%
%   [newexcludelength] returns the new total number of trials now in BadTrials 
%
%
% It works on datasets which already have badtrials marked, by loading up
% the existing ones, concatenating with what you have in your 'exclude' 
% list, removing duplicates, sorting, and resaving. 
% important to load existing files, find/remove duplicates and sort
% numerically before saving update...
%
% YOU NEED YOUR LIST OF BAD TRIALS IN A CELL ARRAY
%
% For safety, set a max-exclude. If the total number of badtrials trying to
% be written to file is bigger than this, then only the existing badtrials
% will be re-written and the stuff you're trying to add will be ignored.




%% First find and load the ClassFile and find it's length. From this we can
%deduce how many badtrials are already listed...
%Dataset = uigetdir('', 'Select CUT dataset')

ClassFile = '/ClassFile.cls'
ClassFileWrite = strcat(Dataset, ClassFile);

fid = fopen(ClassFileWrite, 'r');
C = textscan(fid,'%s', 'delimiter', '\n');

%% First thing to do is update the 'NUMBER OF TRIALS' on Line 22 - should hold total number of BadTrials
replaceline = 22;
%numlines = 25;
totalNoBadTrials = (length(exclude))
newtext = num2str(totalNoBadTrials);
% look at length of current file
filelength = (length(C{1}));
currentexcludeval = (C{1}(22))

% remove 2 blinks if needed here:........
my_temp_cell1 = regexp(C{1}, '\d*', 'match', 'once');
% Convert the cell to numerical array
my_array1 = str2double(my_temp_cell1);
% Replace all NaN's by zeros
my_array1(isnan(my_array1)) = 0;
C2 = my_array1;
C2end = length(C2)
    if C2(C2end) == 0
       C2 = C2(1:(C2end)-1)
    end
C2end = length(C2)
    if C2(C2end) == 0
       C2 = C2(1:(C2end)-1)
    end

% calculate new file length
currnum = currentexcludeval   
newexcludeval = (str2num(currnum{1}) + (str2num(newtext)))

% sometimes it gets confused...this helps, even though it just repeats the
%above.. :/

    if newexcludeval > 0
        newexcludeval = newexcludeval
        else newexcludeval = (currnum) + str2num(newtext)
    end

new2 = num2str(newexcludeval)
fclose(fid)
% open up some more shit
fid = fopen(ClassFileWrite, 'r')
mydata = cell(1, filelength);
    for k = 1:filelength
    mydata{k} = fgetl(fid);
    end
fclose(fid);


% Change the 'existingexclusions' to string
existingexclusions = mydata(25:end);
my_temp_cell = regexp(existingexclusions, '\d*', 'match', 'once');
% Convert the cell to numerical array
my_array = str2double(my_temp_cell);
% Replace all NaN's by zeros
my_array(isnan(my_array)) = 0;
existingexclusions2 = my_array;


% Now concatenate existing with new list... minus1 from the exclude list, sort & remove duplicates
% NB* we -1 from new trial list (those in exclude) because trials start
% from 0 not 1...
newdata = exclude';
newdata = newdata -1;
totalexclude = vertcat(newdata, existingexclusions2');
total = unique(totalexclude);
total = sort(total);
total2 = total(total>0);
newexcludelength = length(total2);


% now, jumping back to updating the total number of trials (couldn't do
% this above as didn't know about duplicates...) update lines 22 and write
% to file:

mydata{22} = num2str(newexcludelength);
fid = fopen(ClassFileWrite,'W');
fprintf(fid, '%s\n', mydata{:});
fclose(fid);


% Format new list of badtrials - add the plus sign
    for s = 1:length(total2);
     newdata2{s}=strcat('+', (num2str(total2(s))));
    end
 
% Get read to print it to file:
% Print that shit with various formatting;:
startline = 25; % where the list of trials starts in a classfile (should be 25)
lineformat ='%s\n'; %the format 'new line', we'll update this in a minute
newformat = '                  %s\n'; %new format....we'll come back to this
newfilelength = 24 + (length(total2))

% re-open up file so we can load in new list of trials
fid = fopen(ClassFileWrite, 'r+');
mydata2 = cell(1, newfilelength);

    for k = 1:newfilelength;
   mydata2{k} = fgetl(fid); %load into mydata2
    end
fclose(fid);
% Put mydata2 (new list) into position in text file
    for n = 25:newfilelength;
    mydata2(n) = newdata2(n-24);
    end
% save it:
fid = fopen(ClassFileWrite, 'w');
fprintf(fid, lineformat, mydata2{:});
fclose(fid);

% Next we need a way to format (indent) list of trials numbers as they're indented 17 spaces from the left.

%try loading up new doc, remove list of trials to new array and save and
%close. then open new doc as 'append' and append trials list (current held
%in array) onto end with appropriate formatting

fid = fopen(ClassFileWrite, 'r+');
%correctdata = cell(1, newfilelength); 
    for t = 1:newfilelength;
   editdata{t} = fgetl(fid); %new correct data is in editdata. just needs reformatting and re saving...
    end
% rename some shit
editdata2 = editdata(1:24); %put top half of doc (un-indented) into editdata2
trials_bkup = editdata; % good idea to backup existing file...
fid = fopen(ClassFileWrite, 'w'); %write top half of file back to file
fprintf(fid, lineformat, editdata2{:});
fclose(fid);

% now create array of trials again, from trials_bkup, and append that onto
% new file with indented formatting

totlength = length(editdata);
trialdata = trials_bkup(25:totlength)'
notrials = length(trialdata);

% This is for safety: if the number of trials being excluded totals greater
% than 50, then will just re-write the RT-excluded data and ignore blink
% detector as shit must've gone cray...

    if notrials > maxexclude;
        trialdata = existingexclusions;
        elseif notrials <= maxexclude;
        trialdata = trialdata;
    end

% BUT, if shit goes cray, if this happens, need to rewrite the total number of exclusions in
% classfile:
shortfilelength = 24
fid = fopen(ClassFileWrite, 'r');
mydata3 = cell(1, shortfilelength);
    for k = 1:shortfilelength;
    mydata3{k} = fgetl(fid);
    end
fclose(fid);
% only if it passes the magical fucking number 50
    if notrials > 50
    mydata3{22} = num2str(currnum);
       elseif notrials < 50;
    mydata3{22} = num2str(notrials);
    end
% write some more shit
fid = fopen(ClassFileWrite,'W');
fprintf(fid, '%s\n', mydata3{:});
fclose(fid);

% Now, re-write the trial list with correct format..(format defines
% earlier)
    for v = 1:length(trialdata);
    fids = fopen(ClassFileWrite,'a');
    fprintf(fids, newformat, trialdata{v});
% close the file 
    fclose(fid);
    end

lastformat2 = '\n\n'
fid = fopen(ClassFileWrite, 'a');
fprintf(fid, lastformat2);
fclose(fid);
    
    
clc
disp(' ')

disp('All your shit has been done for you.');
disp(' ')
disp('Dataset:');
disp(Dataset);
disp(' ')
disp('Total No Exclusions:')
disp(newexcludelength)
end


% for b = 26:length(mydata2)
%     fid = fopen(ClassFileWrite, 'w');
%     fprintf(fid, newformat, mydata2{b});
%     fclose(fid);
% end
% 
% 
% 
% 
% myformat2 = '\n\n'
% fid = fopen(ClassFileWrite, 'a');
% fprintf(fid, myformat2);
% fclose(fid);