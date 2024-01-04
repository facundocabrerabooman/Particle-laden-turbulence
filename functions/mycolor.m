function mycolor_RGB = mycolor(color1,color2,color3)

num12 = 101;num23 = 101;

if ischar(color1)==1
	r1 = double(hex2dec(color1(2:3)))/255;
	g1 = double(hex2dec(color1(4:5)))/255;
	b1 = double(hex2dec(color1(6:7)))/255;
	color1 = [r1, g1, b1];
    
    if nargin>1
        r2 = double(hex2dec(color2(2:3)))/255;
	    g2 = double(hex2dec(color2(4:5)))/255;
	    b2 = double(hex2dec(color2(6:7)))/255;
	    color2 = [r2, g2, b2];

        if nargin>2
            r3 = double(hex2dec(color3(2:3)))/255;
	        g3 = double(hex2dec(color3(4:5)))/255;
	        b3 = double(hex2dec(color3(6:7)))/255;
	        color3 = [r3, g3, b3];
        end
    end
end

if nargin>2
    R_mat=[linspace(color1(1),color2(1),num12),linspace(color2(1),color3(1),num23)];
    G_mat=[linspace(color1(2),color2(2),num12),linspace(color2(2),color3(2),num23)];
    B_mat=[linspace(color1(3),color2(3),num12),linspace(color2(3),color3(3),num23)];
elseif nargin>1
    R_mat=[linspace(color1(1),color2(1),num12)];
    G_mat=[linspace(color1(2),color2(2),num12)];
    B_mat=[linspace(color1(3),color2(3),num12)];
else
    R_mat=r1;
    G_mat=g1;
    B_mat=b1;
end
mycolor_RGB=[R_mat',G_mat',B_mat'];