function [n_max, n, m, timeel] = decodeSPIHT(in, mode, w, h, l)
% decodeSPIHT - SPIHT decoder
%
% parameters:      in - compressed data
% output:          m - image, n_max - max steps, n - last step, timeel - time elapsed 


%-----------  initialization  ----------------

disp(' ');
disp('SPIHT decoder engaged...');
if mode == 'b'
    disp('in BASE mode (full SOTs)');
elseif mode == 'B'    
    if nargin < 5
        disp('fail - in B mode you need to specify w, h, level');
        return;
    else
        size_x = w;
        size_y = h;
        level = l;
    end        
else
    disp('in DEGRADED TREES mode (non-LL root SOTs)');
end

% settings
if mode == 'B'
    n_max = in(1,1);
    ctr = 2;
else
    size_x = in(1,1);
    size_y = in(1,2);
    n_max = in(1,3);
    level = in(1,4);
    ctr = 5;
end

m = zeros(size_y, size_x);

%----------- LIP, LSP, LIS init ----------------
timeel = 0;

% !!!!!!! color level shift enable HERE
bandsizew = 2^(log2(size_x) - level);
bandsizeh = 2^(log2(size_y) - level);

if(bandsizew < 2 || mod(bandsizew, 2) > 0 || bandsizeh < 2 || mod(bandsizeh, 2) > 0 )
    disp('error, bandsize is too small');
end

%-----------   LIP, LSP, LIS init   ----------------
timeel = 0;
LIP = [];
LIS = [];
% for(wCoord j=0; j < bandSizeH_; ++j) {
% 		for(wCoord i=0; i < bandSizeW_; ++i) {
% 			LIP_.push_back(XY(i,j));
% 			if(!(i % 2 == 0 && j % 2 == 0))
% 				LIS_.push_back(XYT(i,j,typeA));
% 		}
% 	}

if mode == 'b';
    % base mode
    for j=1:bandsizeh
        for i=1:bandsizew
            LIP = [LIP; [i j]];
            if ~(mod(i-1,2)==0 && mod(j-1,2) == 0)
                LIS = [LIS; [i j 0]];
            end
        end
    end
else
    % degraded mode
    if bandsizew * 2 > size_x || bandsizeh * 2 > size_y
        disp('error, WT level too small, bandsize*2 is not available');
        return;
    end
    % fill LIP
    for j=1:bandsizeh*2
        for i=1:bandsizew*2
            LIP = [LIP; [i j]];
            if (i > bandsizew || j > bandsizeh)
                LIS = [LIS; [i j 0]];
            end
        end
    end
end

LSP = [];

n = n_max;

disp('DECODER: header successfully processed, initialization finished');

%-----------   decoding   ----------------
n = n_max;
while (ctr <= size(in,2))
    tic;
    ctr_backup=ctr;
    
    %Sorting Pass
    LIPtemp = LIP; temp = 0;
    
    % LIP pass
    for i = 1:size(LIPtemp,1)  
        temp = temp+1;
        
        % bitstream end condition
        if ctr > size(in,2) 
            return
        end
        
        % pixel in LIP is significant
        if in(1,ctr) == 1 
            ctr = ctr + 1;
            
            % bitstream end condition
            if ctr > size(in,2) 
                return
            end
            
            if in(1,ctr) > 0
                %sign +, thr + 1/2 to m
                m(LIPtemp(i,2),LIPtemp(i,1)) = 2^n + 2^(n-1);
            else
                %sign -, -thr - 1/2 to m
                m(LIPtemp(i,2),LIPtemp(i,1)) = -2^n  - 2^(n-1);  
            end
            
            %add pixel to LSP, remove from LIP
            LSP = [LSP; LIPtemp(i,:)];  
            LIP(temp,:) = []; temp = temp - 1;
        end
        
        %not significant -> another pixel
        ctr = ctr + 1;
    end
    
    LIStemp = LIS; temp = 0; i = 1;
    
    while ( i <= size(LIStemp,1))
        temp = temp + 1;
        
        % more?
        if ctr > size(in,2)
            return
        end
        
        % is there significance?
        if in(1,ctr) == 1 
            
            ctr = ctr + 1;
            basex = LIStemp(i,1)-1; basey = LIStemp(i,2)-1;

            if (basex < bandsizew && basey < bandsizeh && mode == 'b')
                if (mod(basey, 2) == 0 && mod(basex, 2) ~= 0)
                    basex = basex + bandsizew - 1;
                elseif (mod(basey, 2) ~= 0 && mod(basex, 2) == 0)
                    basey = basey + bandsizeh - 1;
                else
                    basex = basex + bandsizew - 1; basey = basey + bandsizeh - 1;
                end
            else
                basex = basex * 2; basey = basey * 2;
            end

            for l=1:4
                if l==2
                    basex = basex + 1;
                elseif l==3
                    basex = basex - 1; basey = basey + 1;
                elseif l==4
                    basex = basex + 1;
                end

                %type A processing
                if LIStemp(i,3) == 0

                    % more?
                    if ctr > size(in,2)
                        return
                    end
                    
                    % significant here!
                    if in(1,ctr) == 1
                        
                        ctr = ctr + 1;
                        
                        % read sign
                        if ctr > size(in,2)
                            return
                        end
                        
                        if in(1,ctr) == 1
                            m(basey+1,basex+1) = 2^n + 2^(n-1);  
                        else
                            m(basey+1,basex+1) = -2^n  - 2^(n-1); 
                        end
                        ctr = ctr + 1;
                        LSP = [LSP; basex+1 basey+1];

                    else
                        ctr = ctr + 1;
                        LIP = [LIP; basex+1 basey+1];
                    end

                else
                    % type B processing, make A from them                    
                    LIS = [LIS; basex+1 basey+1 0];
                    LIStemp = [LIStemp;  basex+1  basey+1 0];
                end
            end

            if (LIStemp(i,3) == 0 && basex*2 < size_x && basey*2 < size_y)
                LIS = [LIS;  LIStemp(i,1)  LIStemp(i,2) 1];
                LIStemp = [LIStemp;  LIStemp(i,1)  LIStemp(i,2) 1];              
            end

            % remove [x,y, typ] from LIS
            LIS(temp,:) = []; temp = temp-1;

        else
            % no significance, output -> 0 and no changes
            ctr = ctr + 1;
        end

        i = i+1;
    end

    ctr_backup2 = ctr;
    
    % Refinement Pass
    temp = 1;
    value = m(LSP(temp,2), LSP(temp,1));
    
    % as long as there are items in LSP compliant with the following
    % condition
    while (abs(value) >= 2^(n+1) & (temp <= size(LSP,1)))
        if ctr > size(in,2)
            return
        end

        % add (subtract) according to bit
        value = value + ((-1)^(in(1,ctr) + 1)) * (2^(n-1))*sign(m(LSP(temp,2),LSP(temp,1))); 
        m(LSP(temp,2),LSP(temp,1)) = value;
        ctr = ctr + 1;
        temp = temp + 1;    
        if temp <= size(LSP,1)
            value = m(LSP(temp,2),LSP(temp,1));
        end
    end
    
    time = toc;
    timeel = timeel + time;
    
    disp(['DECODER: STEP ' num2str(n_max-n) ' finished in ' num2str(time) 's. Processed: (S/R): ' num2str(ctr_backup2-ctr_backup) '/' num2str(ctr-ctr_backup2) 'bits. LIS:' num2str(size(LIS,1)) ', LIP:'  num2str(size(LIP,1))  ', LSP:'  num2str(size(LSP,1)) ]);
   
    n = n - 1;
end

   