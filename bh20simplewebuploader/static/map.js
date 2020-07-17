let map = L.map( 'map', {
  center: [37.0902, -95.7129],  // Default to U.S.A
  minZoom: 3,
  zoom: 0
});
L.tileLayer( 'http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
  attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a>',
  subdomains: ['a','b','c']
}).addTo( map );

let markers = L.markerClusterGroup().addTo(map)


function drawMap(){

// initialize the map on the "map" div with a given center and zoom
var mymap = L.map('mapid').setView([51.505, -0.09], 1);

L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
}).addTo(mymap);

fetch(scriptRoot + "api/getCountByGPS")
    .then(response => {
    console.log(response)
      return response.json();
    })
    .then(data => {

   for (var i=0; i<data.length;i++) {
   gps=data[i]["GPS"].split(" ")
    var circle = L.circle([gps[1], gps[0]], {
    color: 'red',
    fillColor: '#f03',
    fillOpacity: 0.5,
    radius: parseInt(data[i]["count"])  //not working for whatever reason
        }).addTo(mymap);
      }

      });
}



/* This function updates the map with markers
 *
*/
function updateMapMarkers() {
    markers.clearLayers(); // remove all markers
    document.getElementById("results").classList.remove("invisible");
    document.getElementById("loader").classList.add("invisible");
    for (let i = 0; i < data.length; i++) {
        let {"count": fastaCount, GPS, LocationLabel: label } = data[i];
        let coordinates = GPS.split(" ");
        if (!(coordinates == null)) {
            let lat, lon;
            [lon, lat] = coordinates.map(parseFloat);
            let point = L.point()
            let marker = L.marker([lat, lon]);
            marker.bindPopup("<b>" + label + "</b><br/>" + "FastaCount: " +fastaCount);
            markers.addLayer(marker)
        }}
    // Reload the map
    map.invalidateSize();
    document.getElementById("map_view").classList.add("invisible");
    document.getElementById("loader").classList.add("invisible");
}
