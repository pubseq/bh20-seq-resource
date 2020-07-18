
var map = L.map( 'mapid', {
    center: [51.505, -0.09],  // Default to U.S.A
    minZoom: 2,
    zoom: 0
});

L.tileLayer( 'http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> | <a href="http://covid19.genenetwork.org/">COVID-19 PubSeq</a>',
    subdomains: ['a','b','c']
}).addTo(map);


function drawMap(){
    var mymap = map;

    fetch(scriptRoot + "api/getCountByGPS")
        .then(response => {
            console.log(response)
            return response.json();
        })
        .then(data => {
            updateMapMarkers(data);

      });
    document.getElementById("map_view").classList.remove("invisible");
    map.invalidateSize();
}




/*
 * This function updates the map with markers
 */

function updateMapMarkers(data) {
    let markers = L.markerClusterGroup({

        iconCreateFunction: function (cluster) {
            var theseMarkers = cluster.getAllChildMarkers(); //// -- this is the array of each marker in the cluster.

            sumCount = 0;
            for (var i = 0; i < theseMarkers.length; i++) {
                // sumCount += theseMarkers[i].feature.properties.Confirmed;
                sumCount += 1;
            }

            var digits = (sumCount + '').length;

            return L.divIcon({
                html: sumCount,
                className: 'cluster digits-' + digits,
                iconSize: null
            });
        },

        pointToLayer: function (feature, latlng) {
            return L.circleMarker(latlng, {
                opacity: 1,
                color: getSColor(feature.properties.Confirmed),
                weight: getSwieght(feature.properties.Confirmed),
                fillColor: getColor(feature.properties.Confirmed),
                fillOpacity: .3,
                radius: getRad(feature.properties.Confirmed),
                pane: 'circlesIIOM'
            });

        },

    });
    for (let i = 0; i < data.length; i++) {
        let {"count": fastaCount, GPS, LocationLabel: label } = data[i];
        let coordinates = GPS.split(" ");
        if (!(coordinates == null)) {
            let lat, lon;
            [lon, lat] = coordinates.map(parseFloat);
            let point = L.point()
            marker = (L.marker([lat, lon]));
            marker.bindPopup("<b>" + label + "</b><br/>" + "SARS-CoV-2<br/>sequences: " +fastaCount);
            markers.addLayer(marker);
        }
    }
    map.addLayer(markers);
}
