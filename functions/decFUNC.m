function dec = decFUNC(X,LoD,HiD,dwtEXTM)

dec = cell(2,2);
permVect = [];
[a_Lo,d_Hi] = wdec1D(X,LoD{1},HiD{1},permVect,dwtEXTM);
permVect = [2,1,3];
[dec{1,1},dec{1,2}] = wdec1D(a_Lo,LoD{2},HiD{2},permVect,dwtEXTM);
[dec{2,1},dec{2,2}] = wdec1D(d_Hi,LoD{2},HiD{2},permVect,dwtEXTM);