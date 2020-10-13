function r = ssimMean(img, x, y, x2, y2)

sum = 0.0;
for j=y:y2
    for i=x:x2
        sum = sum + img(j,i);
    end
end

r = sum / ((x2-x+1) * (y2-y+1));

return;