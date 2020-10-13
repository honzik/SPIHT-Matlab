function value = checkDescendant(i, j, type, m, bsizew, bsizeh, thr, mode)
% checkDescendant - significance search, return 0 - non sign; 1 -
% singificant
%
% parameters:    i,j          - coords of the explored coef
%                type         - A(0) or B(1) type entry
%                m            - tile
%                bsizew       - bandsize w
%                bsizeh       - bandsize h
%                thr          - threshold
%                mode         - mode 'd' or 'b'

if nargin<8, mode='b'; end
value=0;
% due to matlab
x=i-1; y=j-1;
sizew = 2;

if (x<bsizew && y<bsizeh && mode=='b')
    baseX = floor(x / 2.0) * 2;
    baseY = floor(y / 2.0) * 2;
    if mod(y,2) == 0
        baseX = baseX + bsizew;
    else
        if mod(x,2) == 0
            baseY = baseY + bsizeh;
        else
            baseX = baseX + bsizew; baseY = baseY + bsizeh;
        end
    end        
else
    baseX = 2*x; baseY = 2*y;
end

typ = type;

while( baseX < size(m,2) && baseY < size(m,1) )
    if typ == 0
        maxTest = m((baseY+1) : (baseY+sizew) , (baseX+1) : (baseX+sizew));
        if max(max(abs(maxTest))) >= thr
            value=1;
            return;
        end
    else
        typ = 0;
    end
    baseX = 2*baseX; baseY = 2*baseY; sizew = sizew * 2;
end    

    