/*

Draws the map using Leaflet and OpenStreetmap

drawMap() is the main function.

*/
var map = L.map( 'mapid', {
    center: [51.505, -0.09],  // Default to U.S.A
    minZoom: 2,
    zoom: 0
});

L.tileLayer( 'http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> | <a href="http://covid19.genenetwork.org/">COVID-19 PubSeq</a>',
    subdomains: ['a','b','c']
}).addTo(map);

/*
 * When a page gets rendered this function draws the map
 */

function drawMap(){
    var mymap = map;

    // ---- fetch all counts
    fetch(scriptRoot + "api/getCountByGPS")
        .then(response => {
            console.log(response)
            return response.json();
        })
        .then(data => {
            buildMapMarkers(data);

      });
    document.getElementById("map_view").classList.remove("invisible");
    map.invalidateSize();
}


/*
 * Register a marker with special attribute track # sequences
 */

seqMarker = L.Marker.extend({
   options: {
       seqMarkerLocation: "Loc",
       contributors: "countContrib",
       sequences: "countSeq"
   }
});

/*
 * Builds markers on the map. We use cluster groups to allow
 * counts at different zoom levels. This function is called
 * once on page loading. markerClusterGroup just handles it.
 * Note the display is handled in CSS (main.css) as .my-custom-icon*
 */

function buildMapMarkers(data) {
    let markers = L.markerClusterGroup({
        singleMarkerMode: true,
        iconCreateFunction: function (cluster) {
            // ---- add marker
            // array of each marker in the cluster:
            var theseMarkers = cluster.getAllChildMarkers();

            // --- compute zoom level and set style

            sumCount = 0;
            for (var i = 0; i < theseMarkers.length; i++) {
                sumCount += theseMarkers[i].options.sequences;
            }

            if (theseMarkers.length < 2) {
                return L.divIcon({
                    html: sumCount,
                    className: 'my-custom-icon my-custom-icon-0',
                })
            } else {
                var digits = (sumCount + '').length;
                return L.divIcon({
                    html: sumCount,
                    className: 'my-custom-icon my-custom-icon-'+digits,
                });
            }}});
    // ---- Build the marker list
    for (let i = 0; i < data.length; i++) {
        let {"count": fastaCount, GPS, LocationLabel: label } = data[i];
        let countSeq = Number(fastaCount);

        let coordinates = GPS.split(" ");
        if (!(coordinates == null)) {
            let lat, lon;
            [lon, lat] = coordinates.map(parseFloat);
            let point = L.point()
            marker = new seqMarker([lat, lon],markerOptions={title: fastaCount+" sequences",sequences: countSeq});
            marker.bindPopup("<b>" + label + "</b><br/>" + "SARS-CoV-2<br/>sequences: " +fastaCount);
            markers.addLayer(marker);
        }
    }
    map.addLayer(markers);
}
