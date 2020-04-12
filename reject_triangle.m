%%
% returns distance(centriod, circumcenter)/circumcenter_radius.
% consider reject triangle if ratio close to 1
%%
function ratio = reject_triangle(p1, p2, p3)
    point_list = [p1; p2; p3];
    conn_list = [1,2,3];
    tr = triangulation(conn_list, point_list);
    [c,bigR] = circumcenter(tr);
    centriod = (p1+p2+p3)/3;
    smallR = sqrt(dot(centriod-c, centriod-c));
    ratio = smallR/bigR;
end