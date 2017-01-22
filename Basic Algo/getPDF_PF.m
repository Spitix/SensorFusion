function [y, x] = getPDF_PF( xp, dx, xmin, xmax )

    N       =   length( xp ) ; 
    tblx    =   xmin : dx : xmax ;  
    ndata   =   length( tblx ) ; 
    pdf     =   zeros( ndata, 1 ) ; 
    area    =   0.0 ;

    for i = 1 : ndata-1
       
        bin1=   tblx( i ) ; 
        bin2=   tblx(i+1) ; 
        
        for j = 1 : N
            
            if ( bin1 <= xp(j) ) && ( xp(j) < bin2 )                      % Assign pdf
                        
                pdf(i) =   pdf(i) + 1 ;
                        
            end
            
        end
        
        area   =   area + pdf(i) * dx ; 
        
    end
    
    pdf         =   pdf / area ; 
    
    y           =   pdf ; 
    x           =   tblx ; 