function [n_max, n, out, timeel] = encodeSPIHT(m, max_bits, level, mode)
% encodeSPIHT - SPIHT encoder
%
% parameters:        m - source image
%                    max_bits - maximal bit count on the output
%                    level - transform depth
%                    mode - 'd', 'b', 'bb'
% outputs:           out - bitstream,
%                    n_max - max. step count
%                    n - last pass info
%                    timeel - time elapsed

%-----------   initialization  -----------------

disp(' ');
disp('SPIHT enkoder engaged');
if mode == 'b' || mode == 'B'
    disp('in BASE mode (full SOTs)');
elseif mode == 'd'
    disp('in DEGRADED TREES mode (non-LL root SOTs)');
end

bitctr = 0;
size_x = size(m, 2);
size_y = size(m, 1);

out = [];
n_max = floor(log2(max(max(abs(m)))));
Bits_Header = 0;

%-----------   stream header  ----------------
% image size, max_pass and decomp. depth
%                                                                 2*size
%                                                                 level
%                                                                 n_max
if mode=='B'
    out(1,[1]) = [n_max]; bitctr = bitctr + 8; 
    index = 2;
else
    out(1,[1 2 3 4]) = [size_x size_y n_max level]; bitctr = bitctr + 2*16+8+8;
    index = 5;
end

Bits_Header = bitctr;

% compute bandsize

% !!!!!!! color level shift enable HERE
bandsizew = 2^(log2(size_x) - level);
bandsizeh = 2^(log2(size_y) - level);

if(bandsizew < 2 || mod(bandsizew, 2) > 0 || bandsizeh < 2 || mod(bandsizeh, 2) > 0 )
    disp('error, bandsize is too small');
    return
end


%-----------   LIP, LSP, LIS init   ----------------
timeel = 0;
LIP = [];
LIS = [];

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

disp('ENCODER: init finished');
disp(['ENKODER: Max steps = ' num2str(n_max)]);

%-----------   encoding   ----------------
while(n > 0)
    bitctr_backup = bitctr;
    tic;

    % Sorting Pass
    % LIP pass
    LIPtemp = LIP; temp = 0;
    for i = 1:size(LIPtemp,1)
        temp = temp+1;
        if (bitctr + 1) >= max_bits
            % bitstream end condition
            if (bitctr < max_bits)
                out(length(out))=[];
            end
            return
        end

        % significance?
        if abs(m(LIPtemp(i,2),LIPtemp(i,1))) >= 2^n
            % output -> 1
            out(index) = 1; bitctr = bitctr + 1;
            index = index +1;

            % sign
            sgn = m(LIPtemp(i,2),LIPtemp(i,1))>=0;
            out(index) = sgn; bitctr = bitctr + 1;
            index = index +1;

            % into LSP
            LSP = [LSP; LIPtemp(i,:)];

            % LIP remove
            LIP(temp,:) = []; temp = temp - 1;
        else
            % output -> 0 (not significant)
            out(index) = 0; bitctr = bitctr + 1;
            index = index +1;
        end
    end


    
    % LIS pass
    LIStemp = LIS; temp = 0; i = 1;
    while ( i <= size(LIStemp,1))
        temp = temp + 1;
        % LIS is entry type A
        if bitctr >= max_bits
            return
        end

        % returns highest value from the subtree
        signif = checkDescendant(LIStemp(i,1),LIStemp(i,2),LIStemp(i,3),m,bandsizew,bandsizeh,2^n,mode);

        % is there significance?
        if signif==1

            out(index) = 1; bitctr = bitctr + 1;
            index = index +1;
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
                    if (bitctr + 1) >= max_bits
                        if (bitctr < max_bits)
                            out(length(out))=[];
                        end
                        return
                    end

                    if abs(m(basey+1, basex+1)) >= 2^n

                        out(index) = 1; bitctr = bitctr + 1; index = index +1;
                        sgn = m(basey+1, basex+1) >= 0;

                        if (bitctr + 1) >= max_bits
                            if (bitctr < max_bits)
                                out(length(out))=[];
                            end
                            return
                        end

                        out(index) = sgn; bitctr = bitctr + 1;  index = index +1;
                        LSP = [LSP; basex+1 basey+1];

                    else
                        out(index) = 0; bitctr = bitctr + 1; index = index +1;
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
            out(index) = 0; bitctr = bitctr + 1;
            index = index +1;
        end

        i = i+1;


        % type B entry in LIS
        if bitctr >= max_bits
            return
        end

    end
    
    bitctr_backup2 = bitctr;

    % Refinement Pass
    temp = 1;
    value = floor(abs(2^(n_max-n+1)*m(LSP(temp,2),LSP(temp,1))));

    % as long as there are items in LSP compliant with the following
    % condition
    while (value >= 2^(n_max+2) & (temp <= size(LSP,1)))
        if bitctr >= max_bits
            return
        end
        s = bitget(value,n_max+2);
        % save MSB
        out(index) = s; bitctr = bitctr + 1;
        index = index +1;
        temp = temp + 1;
        if temp <= size(LSP,1)
            value = floor(abs(2^(n_max-n+1)*m(LSP(temp,2),LSP(temp,1))));
        end
    end

    time = toc;
    timeel = timeel + time;
    disp(['ENCODER: Step ' num2str(n_max-n) ' finished in ' num2str(time) 's. Sent (S/R): ' num2str(bitctr_backup2-bitctr_backup) '/' num2str(bitctr-bitctr_backup2) 'bits. LIS:' num2str(size(LIS,1)) ', LIP:'  num2str(size(LIP,1))  ', LSP:'  num2str(size(LSP,1)) ]);
    n = n - 1;
end

