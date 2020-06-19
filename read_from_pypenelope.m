clear;
clc;

MAX_TRAG_SEG = 2000;
MAX_PARTICLE = 1000;
%%Load data
fid = fopen('pe-trajectories.dat');
line = fgetl(fid);
trajs = zeros(MAX_TRAG_SEG,4,MAX_PARTICLE);
while line(1,2) == '#'
    line = fgetl(fid);
end
line_cnt = 1;
particle_cnt = 0;
eof = 0;
while ischar(line)
    %fprintf(' Line #%d of text = %s\n', line_cnt, line);
    %fprintf('%s\n',line(1:1));
    if strcmp(line ,'00000000000000000000000000000000000000000000000000000000000000000000000000000000')
        particle_cnt = particle_cnt + 1;
        line_cnt = 1;
        for i = 1 : 7
            if ischar(line)
                line = fgetl(fid);
            else
                eof = 1;
                break;
            end
        end
    end
    if eof == 1
        break;
    end
    aa =  str2num(line);  
    if aa(6) == 1
        for i = 1: 4
            trajs(line_cnt, i ,particle_cnt) = aa(i);
        end
        line_cnt = 1+ line_cnt;
    end
    line = fgetl(fid);
end
save('post_trajs.mat','trajs','particle_cnt')