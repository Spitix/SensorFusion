function [y, x] = getPDF_EKF( mean, P, dx, xmin, xmax )

    tblx    =   xmin : dx : xmax ;  

    pdf     =   ( 1 / sqrt(P) / sqrt(2*pi) ) .*...                  % Assign pdf
                                    exp( - ( tblx - mean ).^2 / 2 / P ) ;  
    
    y           =   pdf ; 
    x           =   tblx ; 


        