
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
 * These functions updates the map with markers
 */

seqMarker = L.Marker.extend({
   options: {
       seqMarkerLocation: "Loc",
       contributors: "countContrib",
       sequences: "countSeq"
   }
});


function updateMapMarkers(data) {
    let markers = L.markerClusterGroup({

        iconCreateFunction: function (cluster) {
            var theseMarkers = cluster.getAllChildMarkers(); //// -- this is the array of each marker in the cluster.

            sumCount = 0;
            for (var i = 0; i < theseMarkers.length; i++) {
                console.log(theseMarkers[i]);
                sumCount += theseMarkers[i].options.sequences;
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
                color: getSColor(10),
                weight: getSwieght(10),
                fillColor: getColor(10),
                fillOpacity: .3,
                radius: getRad(10),
                pane: 'circlesIIOM'
            });

        },

    });
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
