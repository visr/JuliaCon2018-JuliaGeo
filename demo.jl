# the demo was given making use of julia 0.6.4

import ArchGDAL
using PerceptualColourMaps
using GeoInterface
using ColorTypes

GDAL.allregister()

#--- From the lake rough bbox, get image

jsonstring = """{"type":"Polygon","coordinates":[[[-3.01379,54.52648],[-3.00396,54.52648],[-3.00396,54.53042],[-3.01379,54.53042],[-3.01379,54.52648]]]}"""

geom = ArchGDAL.importEPSG(4326) do wgs84
    ArchGDAL.importEPSG(27700) do osgb1936
        ArchGDAL.createcoordtrans(wgs84, osgb1936) do transform
            geom = ArchGDAL.fromJSON(jsonstring)
            ArchGDAL.transform!(geom, transform)
end end end

ArchGDAL.issimple(geom)

geotransform = ArchGDAL.read("lake-district/ortho-20cm.tif") do dataset
    ArchGDAL.getgeotransform(dataset)
end

function rasterindex(centerpoint, geotransform)
    x_min = geotransform[1]
    x_res = geotransform[2]
    y_max = geotransform[4]
    y_res = geotransform[6]
    x, y = centerpoint
    row = Int(fld(y - y_max, y_res)) + 1
    col = Int(fld(x - x_min, x_res)) + 1
    row, col
end

function rasterindex(bbox::GDAL.OGREnvelope, geotransform)
    topleft = (bbox.MinX, bbox.MaxY)
    bottomright = (bbox.MaxX, bbox.MinY)
    r1, c1 = rasterindex(topleft, geotransform)
    r2, c2 = rasterindex(bottomright, geotransform)
    r1:r2, c1:c2
end

bbox = ArchGDAL.getenvelope(geom)
rows, cols = rasterindex(bbox, geotransform)

subortho = ArchGDAL.read("lake-district/ortho-20cm.tif") do dataset
    ArchGDAL.imread(RGB, dataset, (1,2,3), rows, cols)
end

collect(subortho)


#--- Get the three faced peak from the DTM

# indices of peak
rows = 3001:4500
cols = 4001:7000

# read DTM of this area
subdem = ArchGDAL.read("lake-district/dtm-1m.tif") do dataset
    PermutedDimsArray(ArchGDAL.read(dataset, 1, rows, cols), (2,1))
end

# rolling down the hill
function roll!(canvas, dem, i, j)
    count = 1
    while checkbounds(Bool, dem, i-1:i+1, j-1:j+1) && count < 100_000
        surroundings = dem[i-1:i+1, j-1:j+1]
        _, linidx = findmin(surroundings)
        local_i, local_j = ind2sub(surroundings, linidx)
        i += local_i - 2
        j += local_j - 2
        count += 1
        canvas[i-1:i+1, j-1:j+1] = 0.0
    end
end

# hillshade as background map
hillshade = Gray.(relief(collect(subdem), 45.0, 45.0, 1.0))
canvas = copy(hillshade)
roll!(canvas, subdem, 1215, 755); canvas
roll!(canvas, subdem, 1175, 855); canvas
roll!(canvas, subdem, 1240, 800); canvas

#---
