function plotM(m, blocks)

figure(4);

grids{1} = 'ko-';
grids{2} = 'kv:';
grids{3} = 'k*-.';
grids{4} = 'k+--';
grids{5} = 'm+-';
grids{6} = 'y+-';
grids{7} = 'k+-';

len = size(blocks,2);

if(len > 7)
    disp('plotM message: Too many plots, max 7 simultaneous allowed, truncating...');
	len=7;
    blocks=blocks(1,1:7);
end


switch(len)
    case 1
        plot(m{blocks(1,1)}(2,:), m{blocks(1,1)}(1,:), grids{1});
    case 2
        plot(m{blocks(1,1)}(2,:), m{blocks(1,1)}(1,:), grids{1}, m{blocks(1,2)}(2,:), m{blocks(1,2)}(1,:), grids{2});
    case 3
        plot(m{blocks(1,1)}(2,:), m{blocks(1,1)}(1,:), grids{1}, m{blocks(1,2)}(2,:), m{blocks(1,2)}(1,:), grids{2}, m{blocks(1,3)}(2,:), m{blocks(1,3)}(1,:), grids{3});
    case 4
        plot(m{blocks(1,1)}(2,:), m{blocks(1,1)}(1,:), grids{1}, m{blocks(1,2)}(2,:), m{blocks(1,2)}(1,:), grids{2}, m{blocks(1,3)}(2,:), m{blocks(1,3)}(1,:), grids{3}, m{blocks(1,4)}(2,:), m{blocks(1,4)}(1,:), grids{4});
    case 5
        plot(m{blocks(1,1)}(2,:), m{blocks(1,1)}(1,:), grids{1}, m{blocks(1,2)}(2,:), m{blocks(1,2)}(1,:), grids{2}, m{blocks(1,3)}(2,:), m{blocks(1,3)}(1,:), grids{3}, m{blocks(1,4)}(2,:), m{blocks(1,4)}(1,:), grids{4}, m{blocks(1,5)}(2,:), m{blocks(1,5)}(1,:), grids{5} );
    case 6
        plot(m{blocks(1,1)}(2,:), m{blocks(1,1)}(1,:), grids{1}, m{blocks(1,2)}(2,:), m{blocks(1,2)}(1,:), grids{2}, m{blocks(1,3)}(2,:), m{blocks(1,3)}(1,:), grids{3}, m{blocks(1,4)}(2,:), m{blocks(1,4)}(1,:), grids{4}, m{blocks(1,5)}(2,:), m{blocks(1,5)}(1,:), grids{5}, m{blocks(1,6)}(2,:), m{blocks(1,6)}(1,:), grids{6} );
    case 7
        plot(m{blocks(1,1)}(2,:), m{blocks(1,1)}(1,:), grids{1}, m{blocks(1,2)}(2,:), m{blocks(1,2)}(1,:), grids{2}, m{blocks(1,3)}(2,:), m{blocks(1,3)}(1,:), grids{3}, m{blocks(1,4)}(2,:), m{blocks(1,4)}(1,:), grids{4}, m{blocks(1,5)}(2,:), m{blocks(1,5)}(1,:), grids{5}, m{blocks(1,6)}(2,:), m{blocks(1,6)}(1,:), grids{6}, m{blocks(1,7)}(2,:), m{blocks(1,7)}(1,:), grids{7});
end




