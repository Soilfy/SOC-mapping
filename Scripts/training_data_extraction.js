// Load the CSV file
var bihar_points = ee.FeatureCollection('projects/ee-najah/assets/cimmyt_bihar1');
Map.centerObject(bihar_points);
Map.addLayer(bihar_points, {}, "{Points");

print(bihar_points.limit(100));

//get bihar shapefile


var bihar = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', 'Bihar'));

Map.addLayer(bihar, {}, "Bihar shapefile");

// get sentinel

//date range

var startDate = '2019-11-01';
var endDate = '2019-11-30';

//cloud masking


/**
* Function to mask clouds using the Sentinel-2 QA band
* @param {ee.Image} image Sentinel-2 image
* @return {ee.Image} cloud masked Sentinel-2 image
*/
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000);
}

//filter sentinel


var s2 = ee.ImageCollection('COPERNICUS/S2_SR')
          .filterDate(startDate, endDate)
          .filterBounds(bihar)
          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',10))
          .map(maskS2clouds);
var s2_bihar = s2.median().clip(bihar);
          
          
//print(s2.size());
// Select the bands of interest
var bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];

// // Clip the images to the Bihar state boundary
// var sentinel_bihar = s2.select(bands)
//                             .map(function(image) {
//                               return image.clip(bihar);
//                             });
                            
                            
//viz

// Define the visualization parameters
var visParams = {
  bands: ['B4', 'B3', 'B2'],
  min: 0,
  max: 0.3,
  gamma: 1.4
};

// Add the imagery to the map
Map.addLayer(s2_bihar, visParams, 'Sentinel-2 Bihar');
// elevation and slope


// Load the SRTM dataset
var srtm = ee.Image('USGS/SRTMGL1_003');

// Clip the SRTM dataset to the Bihar shapefile
var srtm_bihar = srtm.clip(bihar);

// Compute the elevation and slope from the SRTM dataset
var elevation = srtm_bihar.select('elevation');
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

s2_bihar = s2_bihar.addBands([elevation_10m, slope_10m]);

// add lst


// add lst

// Load the land surface temperature dataset
var dataset = ee.ImageCollection('MODIS/061/MOD11A1').filter(ee.Filter.date(startDate, endDate));
var landSurfaceTemperature = dataset.select('LST_Day_1km').mean().clip(bihar);

// Downsample the land surface temperature band to 10 meters
var lst_10m = landSurfaceTemperature.resample('bilinear').reproject({
  'crs': landSurfaceTemperature.projection(),
  'scale': 10
});

// Add the land surface temperature band to the image
var s2_bihar = s2_bihar.addBands(lst_10m);


// add to the lsit


bands = bands.concat(['elevation', 'slope', 'LST_Day_1km']);


//create indices

// Compute the NDVI index and add it as a band to the image.
var ndvi = s2_bihar.expression('(nir - red) / (nir + red)', {
  'nir': s2_bihar.select('B8'), // NIR (near-infrared) band
  'red': s2_bihar.select('B4')  // Red band
}).rename('ndvi');

// Compute the NDWI index and add it as a band to the image.
var ndwi = s2_bihar.expression('(nir - green) / (nir + green)', {
  'nir': s2_bihar.select('B8'), // NIR band
  'green': s2_bihar.select('B3') // Green band
}).rename('ndwi');

// Compute the NDMI index and add it as a band to the image.
var ndmi = s2_bihar.expression('(nir - swir1) / (nir + swir1)', {
  'nir': s2_bihar.select('B8'),   // NIR band
  'swir1': s2_bihar.select('B11') // SWIR1 band
}).rename('ndmi');

// Compute the NBR index and add it as a band to the image.
var nbr = s2_bihar.expression('(nir - swir2) / (nir + swir2)', {
  'nir': s2_bihar.select('B8'),   // NIR band
  'swir2': s2_bihar.select('B12') // SWIR2 band
}).rename('nbr');

// Compute the BAI index and add it as a band to the image.
var bai = s2_bihar.expression('1.0 / ((0.1 - red) ** 2 + (0.06 - nir) ** 2)', {
  'red': s2_bihar.select('B4'),  // Red band
  'nir': s2_bihar.select('B8')   // NIR band
}).rename('bai');

// Compute the BSI index and add it as a band to the image.
var bsi = s2_bihar.expression('(swir2 + red) / (nir + swir2)', {
  'swir2': s2_bihar.select('B12'), // SWIR2 band
  'red': s2_bihar.select('B4'),    // Red band
  'nir': s2_bihar.select('B8')     // NIR band
}).rename('bsi');


// Compute the EVI index and add it as a band to the image.
var evi = s2_bihar.expression('2.5 * ((nir - red) / (nir + 6 * red - 7.5 * blue + 1))', {
  'nir': s2_bihar.select('B8'), // NIR band
  'red': s2_bihar.select('B4'), // Red band
  'blue': s2_bihar.select('B2') // Blue band
}).rename('evi');

// Compute the SAVI index and add it as a band to the image.
var savi = s2_bihar.expression('1.5 * ((nir - red) / (nir + red + 0.5))', {
  'nir': s2_bihar.select('B8'), // NIR band
  'red': s2_bihar.select('B4')  // Red band
}).rename('savi');

// Compute the LAI index and add it as a band to the image.
var lai = s2_bihar.expression('3.618 * (nir - red) / (nir + red + 0.684)', {
  'nir': s2_bihar.select('B8'), // NIR band
  'red': s2_bihar.select('B4')  // Red band
}).rename('lai');

// Compute the Albedo and add it as a band to the image.
var albedo = s2_bihar.expression('(1.0 - ((nir + red) - (swir1 + swir2)) / (nir + red + swir1 + swir2))', {
  'nir': s2_bihar.select('B8'),   // NIR band
  'red': s2_bihar.select('B4'),   // Red band
  'swir1': s2_bihar.select('B11'), // SWIR1 band
  'swir2': s2_bihar.select('B12')  // SWIR2 band
}).rename('albedo');

// Update the 'bands' list to include these additional indices
bands = bands.concat(['ndvi', 'ndwi', 'ndmi', 'nbr', 'bai', 'bsi', 'savi', 'evi', 'lai']);

// Add the computed indices to the image
s2_bihar = s2_bihar.addBands([ndvi, ndwi, ndmi, nbr, bai, bsi, savi, evi, albedo, lai]);

// trainig data creation


var training_data = s2_bihar.select(bands).sampleRegions({
  'collection': bihar_points,
  'properties': ['OC'],
  'scale': 10
});

print(training_data.size());

print(training_data.limit(100));

// split the data to test and train

var training_data = training_data.randomColumn();

//export it 


Export.table.toDrive({
  collection: training_data,
  description: 'spolify_bihar_training',
  fileFormat: 'csv',
  folder: 'ee'
});

