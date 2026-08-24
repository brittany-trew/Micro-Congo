// ========== USER CONFIGURATION ==========
var siteLabel = 'Imbalanga';  // << change this once
var utmCRS = 'EPSG:32633';       // << set once for the region
var exportFolder = 'MODIS_LAI_' + siteLabel;
var assetPath = 'projects/ee-brittanytrew01/assets/' + siteLabel + '_bbox';
// =========================================

// Load bounding box as FeatureCollection
var bbox = ee.FeatureCollection(assetPath);
var extentGeometry = bbox.geometry();

// Load MODIS LAI collection
var modisLAI = ee.ImageCollection("MODIS/061/MCD15A3H")
  .select('Lai')
  .filterBounds(extentGeometry);

// Load bounding box asset
var bbox = ee.FeatureCollection(assetPath);
Map.centerObject(bbox, 8);
Map.addLayer(bbox, {}, 'Site Bounding Box');


// Define years
var startYear = 2004;
var endYear = 2024;
var years = ee.List.sequence(startYear, endYear);

// Processing function for one year
var exportYearlyLAI = function(year) {
  year = ee.Number(year);
  var yearStart = ee.Date.fromYMD(year, 1, 1);
  var yearEnd = yearStart.advance(1, 'year');

  var yearlyCollection = modisLAI.filterDate(yearStart, yearEnd);

  var monthlyMeans = ee.List.sequence(1, 12).map(function(m) {
    var month = ee.Number(m);
    var monthStart = ee.Date.fromYMD(year, month, 1);
    var monthEnd = monthStart.advance(1, 'month');

    return yearlyCollection.filterDate(monthStart, monthEnd).mean()
      .set('month', month);
  });

  var yearlyImage = ee.ImageCollection.fromImages(monthlyMeans)
    .toBands()
    .clip(extentGeometry)
    .reproject({
      crs: utmCRS,
      scale: 500
    })
    .set('system:time_start', yearStart.millis());

  Export.image.toDrive({
    image: yearlyImage,
    description: 'MODIS_LAI_Year_' + year.getInfo() + '_' + siteLabel,
    folder: exportFolder,
    fileNamePrefix: 'MODIS_LAI_' + siteLabel + '_' + year.getInfo(),
    region: extentGeometry,
    scale: 500,
    crs: utmCRS,
    maxPixels: 1e13
  });
};

// Loop over all years
years.getInfo().forEach(exportYearlyLAI);
