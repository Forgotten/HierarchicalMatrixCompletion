function prob = sampling_ratio(m)
 
switch m  
    case 4
        prob = 1;
    case 8
        prob = 1;
    case 16 
        prob = 0.90;   
    case 32
        prob = 0.60;
    case 64
        prob = 0.25;
    case 128 
        prob = 0.15;

end