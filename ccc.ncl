load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
level = (/850,700,500/)
limit = (/108,126.5,23,42.3/)  ;west,east,south,north
path1     = "."
path2     = "high/plot/"
path3     = "2008"
spl_name  = tochar(path3)
path4     = spl_name(0:1)+"_"+spl_name(2:3)

filename1 = path1+"before"+"/wrfout_d01_2016-01-"+path4+"0000"
filename2 = path1+"after"+"/wrfout_d01_2016-01-"+path4+"0000"
a = addfile(filename1,"r")
b = addfile(filename2,"r")

lat_before = a->XLAT(0,:,:)
lon_before = a->XLONG(0,:,:)
lat_after  = b->XLAT(0,:,:)
lon_after  = b->XLONG(0,:,:)
z_before   = wrf_user_getvar(a, "z",0)
tc_before  = wrf_user_getvar(a,"tc",0)
u_before   = wrf_user_getvar(a,"ua",0)
v_before   = wrf_user_getvar(a,"va",0)
p_before   = wrf_user_getvar(a, "pressure",0)

z_after    = wrf_user_getvar(b, "z",0)
tc_after   = wrf_user_getvar(b,"tc",0)
u_after    = wrf_user_getvar(b,"ua",0)
v_after    = wrf_user_getvar(b,"va",0)
p_after    = wrf_user_getvar(b, "pressure",0)

kl        = dimsizes(level)
do pl = 0,kl-1

    pressure = level(pl)
    tc_plane_before = wrf_user_intrp3d(tc_before,p_before,"h",pressure,0.,False)
    z_plane_before  = wrf_user_intrp3d( z_before,p_before,"h",pressure,0.,False)/10.
    u_plane_before  = wrf_user_intrp3d( u_before,p_before,"h",pressure,0.,False)
    v_plane_before  = wrf_user_intrp3d( v_before,p_before,"h",pressure,0.,False)

    tc_plane_after  = wrf_user_intrp3d(tc_after,p_after,"h",pressure,0.,False)
    z_plane_after   = wrf_user_intrp3d( z_after,p_after,"h",pressure,0.,False)/10.
    u_plane_after   = wrf_user_intrp3d( u_after,p_after,"h",pressure,0.,False)
    v_plane_after   = wrf_user_intrp3d( v_after,p_after,"h",pressure,0.,False)

    filename3 = path2+tostring(pressure)+"/1601"+path3+".000"

    c = readAsciiTable(filename3,10,"float",2)
    lon        = c(:,1)
    lat        = c(:,2)
    z          = c(:,5)
    tc         = c(:,6)
    direc      = c(:,8)
    spd        = c(:,9)
    uv         = wind_component(spd,direc,1)
    u          = uv[0]
    v          = uv[1]



    aaa = NewList("lifo")
    number = dimsizes(lon)
    do i = 0,number-1
        if lon(i) .ge. limit(0) .and. lon(i) .le. limit(1) .and.lat(i) .ge. limit(2) .and. lat(i) .le. limit(3)
            ListAppend(aaa,i)
        end if
    end do
    in_number = ListCount(aaa)
    if in_number .ne. 0
        in_index = new(in_number,"integer")
        j = 0
        do i = 0,number-1
            if lon(i) .ge. limit(0) .and. lon(i) .le. limit(1) .and.lat(i) .ge. limit(2) .and. lat(i) .le. limit(3)
                in_index(j) = i
                j = j+1
            end if
        end do
        z_before_station         = rcm2points_Wrap(lat_before,lon_before,z_plane_before,lat(in_index),lon(in_index),0)
        z_after_station          = rcm2points_Wrap(lat_before,lon_before,z_plane_after,lat(in_index),lon(in_index),0)

        t_before_station         = rcm2points_Wrap(lat_before,lon_before,tc_plane_before,lat(in_index),lon(in_index),0)
        t_after_station          = rcm2points_Wrap(lat_before,lon_before,tc_plane_after,lat(in_index),lon(in_index),0)

        u_before_station         = rcm2points_Wrap(lat_before,lon_before,u_plane_before,lat(in_index),lon(in_index),0)
        u_after_station          = rcm2points_Wrap(lat_before,lon_before,u_plane_after,lat(in_index),lon(in_index),0)

        v_before_station         = rcm2points_Wrap(lat_before,lon_before,v_plane_before,lat(in_index),lon(in_index),0)
        v_after_station          = rcm2points_Wrap(lat_before,lon_before,v_plane_after,lat(in_index),lon(in_index),0)


        pearson_z_before_station = escorc(z_before_station,z(in_index))
        pearson_z_after_station  = escorc(z_after_station, z(in_index))

        pearson_t_before_station = escorc(t_before_station,tc(in_index))
        pearson_t_after_station  = escorc(t_after_station, tc(in_index))

        pearson_u_before_station = escorc(u_before_station,u(in_index))
        pearson_u_after_station  = escorc(u_after_station, u(in_index))

        pearson_v_before_station = escorc(v_before_station,v(in_index))
        pearson_v_after_station  = escorc(v_after_station, v(in_index))

        print(pearson_z_before_station)
        print(pearson_z_after_station)
        if pearson_z_before_station .gt. pearson_z_after_station
            print("z before is greater")
        else if pearson_z_before_station .lt. pearson_z_after_station
            print("z after is greater")
        else
            print("zs are the same")
        end if
        end if
        print("-------------  "+tostring(pressure)+" z  -------------")

        print(pearson_t_before_station)
        print(pearson_t_after_station)
        if pearson_t_before_station .gt. pearson_t_after_station
            print("t before is greater")
        else if pearson_t_before_station .lt. pearson_t_after_station
            print("t after is greater")
        else
            print("ts are the same")
        end if
        end if
        ;print(pearson_u_before_station)
        ;print(pearson_u_after_station)
        ;print(pearson_v_before_station)
        ;print(pearson_v_after_station)
        print("=============  "+tostring(pressure)+" t  =============")
    else
        print("no value")
    end if

end do
end
