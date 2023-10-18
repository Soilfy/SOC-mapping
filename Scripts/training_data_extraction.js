// Load the CSV file
// Replace the point data (biaf,cimmyt, pnp, wrms) with your data
//  Replace roi with your data 
//var roi = `your_location`


// get sentinel data

var s2_harmonized = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterDate('2020-05-01', '2020-05-31') // *****Check date*****
                  .filterBounds(roi)
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20));
var visualization = {
  min: 0.0,
  max: 3000,
  bands: ['B4', 'B3', 'B2'],
};
var sat_img = s2_harmonized
              .mosaic()
              .clip(roi)
              .select(['B2','B3','B4','B5','B6','B7','B8','B11','B12']);
Map.addLayer(sat_img, visualization, 'RGB');

var imgVV = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filterDate('2020-05-01', '2020-05-31') // *****Check date*****
        .select('VV')
        .map(function(image) {
          var edge = image.lt(-30.0);
          var maskedImage = image.mask().and(edge.not());
          return image.updateMask(maskedImage);
        });

var imgVH = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filterDate('2020-05-01', '2020-05-31') // *****Check date*****
        .select('VH')
        .map(function(image) {
          var edge = image.lt(-30.0);
          var maskedImage = image.mask().and(edge.not());
          return image.updateMask(maskedImage);
        });

var descVV = imgVV.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')).mosaic().clip(roi);
var descVH = imgVH.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')).mosaic().clip(roi);
Map.addLayer(descVV,{min:-15,max:3},'descVV');
Map.addLayer(descVH,{min:-15,max:5},'descVH');

var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_BY_HOUR')
                .filter(ee.Filter.date('2020-05-01', '2020-05-31')); // *****Check date*****

var tempimg = era5.median().clip(roi).select(['skin_temperature','total_precipitation',]);

//create indices

// Compute the NDVI index and add it as a band to the image.
var ndvi = sat_img.expression('(nir - red) / (nir + red)', {
  'nir': sat_img.select('B8'), // NIR (near-infrared) band
  'red': sat_img.select('B4')  // Red band
}).rename('ndvi');

// Compute the NDWI index and add it as a band to the image.
var ndwi = sat_img.expression('(nir - green) / (nir + green)', {
  'nir': sat_img.select('B8'), // NIR band
  'green': sat_img.select('B3') // Green band
}).rename('ndwi');

// Compute the NDMI index and add it as a band to the image.
var ndmi = sat_img.expression('(nir - swir1) / (nir + swir1)', {
  'nir': sat_img.select('B8'),   // NIR band
  'swir1': sat_img.select('B11') // SWIR1 band
}).rename('ndmi');

// Compute the NBR index and add it as a band to the image.
var nbr = sat_img.expression('(nir - swir2) / (nir + swir2)', {
  'nir': sat_img.select('B8'),   // NIR band
  'swir2': sat_img.select('B12') // SWIR2 band
}).rename('nbr');

// Compute the BAI index and add it as a band to the image.
var bai = sat_img.expression('1.0 / ((0.1 - red) ** 2 + (0.06 - nir) ** 2)', {
  'red': sat_img.select('B4'),  // Red band
  'nir': sat_img.select('B8')   // NIR band
}).rename('bai');

// Compute the BSI index and add it as a band to the image.
var bsi = sat_img.expression('(swir2 + red) / (nir + swir2)', {
  'swir2': sat_img.select('B12'), // SWIR2 band
  'red': sat_img.select('B4'),    // Red band
  'nir': sat_img.select('B8')     // NIR band
}).rename('bsi');


// Compute the EVI index and add it as a band to the image.
var evi = sat_img.expression('2.5 * ((nir - red) / (nir + 6 * red - 7.5 * blue + 1))', {
  'nir': sat_img.select('B8'), // NIR band
  'red': sat_img.select('B4'), // Red band
  'blue': sat_img.select('B2') // Blue band
}).rename('evi');

// Compute the SAVI index and add it as a band to the image.
var savi = sat_img.expression('1.5 * ((nir - red) / (nir + red + 0.5))', {
  'nir': sat_img.select('B8'), // NIR band
  'red': sat_img.select('B4')  // Red band
}).rename('savi');

// Compute the LAI index and add it as a band to the image.
var lai = sat_img.expression('3.618 * (nir - red) / (nir + red + 0.684)', {
  'nir': sat_img.select('B8'), // NIR band
  'red': sat_img.select('B4')  // Red band
}).rename('lai');

// Compute the Albedo and add it as a band to the image.
var albedo = sat_img.expression('(1.0 - ((nir + red) - (swir1 + swir2)) / (nir + red + swir1 + swir2))', {
  'nir': sat_img.select('B8'),   // NIR band
  'red': sat_img.select('B4'),   // Red band
  'swir1': sat_img.select('B11'), // SWIR1 band
  'swir2': sat_img.select('B12')  // SWIR2 band
}).rename('albedo');

// Load the SRTM dataset
var srtm = ee.Image('USGS/SRTMGL1_003');

// Clip the SRTM dataset to the Bihar shapefile
var srtm_roi = srtm.clip(roi);

// Compute the elevation and slope from the SRTM dataset
var elevation = srtm_roi.select('elevation');
var slope = ee.Terrain.slope(elevation);

// Downsample the elevation and slope bands to 10 meters
var elevation_10m = elevation.resample('bilinear').reproject({
  crs: elevation.projection(),
  scale: 10
});
var slope_10m = slope.resample('bilinear').reproject({
  crs: slope.projection(),
  scale: 10
});

sat_img = sat_img.addBands([elevation_10m, slope_10m]);

// Update the 'bands' list to include these additional indices
var img = sat_img.addBands([tempimg,descVV,descVH,ndvi,ndwi,ndmi,nbr,bai,bsi,savi,evi,lai]);
print(img);
Map.addLayer(img,{}, 'Concat');
// bands = bands.concat(['ndvi', 'ndwi', 'ndmi', 'nbr', 'bai', 'bsi', 'savi', 'evi', 'lai']);
// ---------------------------------------------------------------------------------------
// -------------------------------- Vector layer -----------------------------------------
// ---------------------------------------------------------------------------------------

// trainig data creation

var full_csv = pnp.merge(biaf).merge(wrms); // They are all done in May 2022
print('Size of csv points:',full_csv.size());
Map.addLayer(full_csv);
Map.centerObject(full_csv,8);

var training_data = img.sampleRegions({
  'collection': full_csv,
  'properties': ['Id','OC'],
  'scale': 10
});

print('Size of sampled points:',training_data.size());

print(training_data.limit(100));

Export.table.toDrive({
  collection: training_data,
  description: 'biaf_wrms_pnp__sample',
  fileFormat: 'csv',
  folder: 'ee'
});
