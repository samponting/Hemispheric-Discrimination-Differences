function MB = RGBtoMB(RGB)
    %Conversion from RGB image to MB coordinate
    %RGB = imread('/Users/samuelponting/Downloads/lightprobe/HS_En1.png');
    load Matrix_LMS2RGB
    
    Lw = 0.689903;
    Mw = 0.348322;
    Sw = 0.0371597/0.0192;

    if length(size(RGB)) == 3
        RGB_v = reshape(RGB,size(RGB,1)*size(RGB,2),size(RGB,3));
    elseif length(size(RGB)) == 2
        RGB_v = RGB;
    end
    
    LMS_v = (Matrix_LMS2RGB*double(RGB_v'))';
        
    MB_v = zeros(size(RGB_v));
    
    MB_v(:,3) = Lw*LMS_v(:,1)+Mw*LMS_v(:,2);
    MB_v(:,1) = Lw*LMS_v(:,1)./MB_v(:,3);
    MB_v(:,2) = Sw*LMS_v(:,3)./MB_v(:,3);
    MB_v(:,3) = MB_v(:,3)/max(MB_v(:,3));
    MB_v(MB_v(:,3)==0,1:2) = 0;
    MB = reshape(MB_v,size(RGB,1),size(RGB,2),size(RGB,3));
    
end