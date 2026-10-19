// ========== USER CONFIGURATION ==========
var siteLabel = 'Imbalanga';  // << change this once
var utmCRS = 'EPSG:32736';       // << set once for the region
//var exportFolder = 'MODIS_LAI_' + siteLabel;
var assetPath = 'projects/ee-username/assets/' + siteLabel + '_bbox';
// =========================================

// Function to mask clouds and shadows using the SCL band
function maskS2sr(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(3) // cloud shadow
              .and(scl.neq(8)) // cloud medium probability
              .and(scl.neq(9)) // cloud high probability
              .and(scl.neq(10)) // thin cirrus
              .and(scl.neq(11)) // snow or ice
              .and(scl.neq(0)) // no data
              .and(scl.neq(1)); // saturated or defective
  return image.updateMask(mask);
}

// Load bounding box
var bbox = ee.FeatureCollection(assetPath);

// Load and mask Sentinel-2 SR images
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(bbox)
  .filterDate('2021-01-01', '2021-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 50)) // still useful for speed
  .map(maskS2sr);

// RGB composite (Red, Green, Blue)
var rgb = s2.select(['B4', 'B3', 'B2'])
  .median()  // or use .qualityMosaic('NDVI') if relevant
  .clip(bbox);

// CIR composite (NIR, Red, Green)
var cir = s2.select(['B8', 'B4', 'B3'])
  .median()
  .clip(bbox);

// Export RGB
Export.image.toDrive({
  image: rgb,
  description: siteLabel+'_RGB_Sentinel2_2021_masked',
  folder: 'GEE_exports',
  fileNamePrefix: siteLabel+'_Sentinel2_RGB_2021_masked',
  region: bbox.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export CIR
Export.image.toDrive({
  image: cir,
  description: siteLabel+'_CIR_Sentinel2_2021_masked',
  folder: 'GEE_exports',
  fileNamePrefix: siteLabel+'_Sentinel2_CIR_2021_masked',
  region: bbox.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
