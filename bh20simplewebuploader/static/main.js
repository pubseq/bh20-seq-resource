/*
 * Menu and navigation
 */

/* Small helpers */
function cell(text) {
    html = ""
    html += '  <div class="rTableCell">';
    html += text;
    html += '  </div>';
    return html
}

/* Convert a list of table items to an HTML DIV table */
function toDIVTable(rows) {
    if (rows.length == 0)
        html = "no results";
    else {
        html = '<div class="rTable">';
        rows.forEach(row => {
            id = row['id'];
            info = row['info'];
            html += '<div class="rTableRow">';
            html += cell('<a href="'+info+'">'+id+'</a>');
            html += cell('<a href="'+row['collection']+'">Collection</a>');
            html += cell('<a href="'+row['fasta']+'">FASTA</a>');
            html += cell('<a href="/api/ebi/sample-'+id+'.xml">EBI/ENA export XML</a>');
            html += '</div>';
        });
        html += '</div>';
    }
    console.log(html);
    document.getElementById("table").innerHTML = html;
}

/* Toggle between adding and removing the "responsive" class to topnav
 * when the user clicks on the icon */
function myFunction() {
    var x = document.getElementById("myTopnav");
    if (x.className === "topnav") {
        x.className += " responsive";
    } else {
        x.className = "topnav";
    }
}

function fetchAPI(apiEndPoint,injectHTML=NULL) {
    fetch(scriptRoot + apiEndPoint)
        .then(response => {
            console.log("* response",response);
            return response.json();
        })
        .then(data => {
            console.log("* data",data);
            if (injectHTML)
                injectHTML(data);
        });
}

// Copy from function above but now output HTML table instead of plain json
function fetchHTMLTable(apiEndPoint) {
    fetch(scriptRoot + apiEndPoint)
        .then(response => {
            return response.json();
        })
        .then(data => {
            console.log(data)
            htmlString="<table>"
            for (var i=0; i<data.length;i++) {
                let url = data[i]["key"];
                let label = data[i]["label"];
                htmlString=htmlString+"<tr><td><a href=\""+url+"\">"+label+"</a></td><td>"+data[i]["count"]+"<td></tr>"
            }
            htmlString=htmlString+"</table>"

            document.getElementById("table").innerHTML = htmlString
        });
}

/* Fetch record info using a 'global search'. Returns for example

[
  {
    "id": "MT326090.1",
    "seq": "http://collections.lugli.arvadosapi.com/c=10eaef75e0b875f81aa1f411c75370cf+126/sequence.fasta"
  },
  {
    "id": "MT326090.1",
    "seq": "http://collections.lugli.arvadosapi.com/c=5a4c815f3e076ad7760a91864c39dd07+126/sequence.fasta"
  }
]

 */
let searchGlobal = (toHTML) => {
    let m =  document.getElementById('search-input').value;
    fetchAPI(scriptRoot + "/api/search?s=" + encodeURIComponent(m), toHTML);
}

// Same as above, but generates div table
let searchGlobaltoDIVTable = () => {
    searchGlobal(toDIVTable);
}

let searchSeq = () => {
  let m =  document.getElementById('search-input').value;
  fetchAPI(scriptRoot + "/api/getDetailsForSeq?seq=" + encodeURIComponent(m));
}

// Get count from Arvados
let fetchCount = () => {
  fetchAPI("/api/getCount");
}

// Get count from Virtuoso
let fetchCountDB = () => {
  fetchAPI("/api/getCountDB");
}

let fetchSEQCountBySpecimen = () => {
  fetchHTMLTable("/api/getSEQCountbySpecimenSource");
}

let fetchSEQCountByLocation = () => {
  fetchHTMLTable("/api/getSEQCountbyLocation");
}

let fetchSEQCountByTech = () => {
  fetchHTMLTable("/api/getSEQCountbytech");
}

let fetchAllaccessions = () => {
  fetchHTMLTable("/api/getAllaccessions");
};

let fetchMap = () => {
    fetchAPI("/api/getCountByGPS");
    updateMapMarkers();
};

let fetchSEQByLocation = () => {
  console.log("Missing - set parameter for request, retrieve data")
};

let fetchSEQCountbyContinent = () => {
  fetchHTMLTable("/api/getSEQCountbyContinent");
}


/*
 * Make sure that only one of the manual metadata entry and metadata upload
 * form components is *actually* a child of the form element in the DOM.
 *
 * Because both make use of the "required" attribute, we can't get away with
 * just hiding the one we don't want the user to fill in. The hidden part will
 * still have possibly empty required fields and (some) browsers will
 * blocksubmission because of it. Moreover, the data (including file uploads)
 * from the hidden elements will still be sent to the server, which the user
 * may not expect.
 */


function setUploadMode() {
  // Make the upload form the one in use.
  uploadFormSpot.appendChild(uploadForm)
  // Remove the upload form from the DOM so its required-ness does not block submission.
  fillFormSpot.removeChild(fillForm)
}

function setFillMode() {
  // Make the fillable form the one in use
  uploadFormSpot.removeChild(uploadForm)
  // Remove the fillable form from the DOM so its required-ness does not block submission.
  fillFormSpot.appendChild(fillForm)
}

function setMode() {
  // Pick mode based on radio
  if (document.getElementById('metadata_upload').checked) {
    setUploadMode()
 } else {
    setFillMode()
 }
}

/*
 * Machinery for variable-length lists of input items.
 */

// Start in mode appropriate to selected form item.
// It is important that we run this code when the page starts! The browser may
// have set the radio button to whatever the state was on last page load,
// instead of the default state, without raising an event, and we have to
// handle that.

/**
 * Add another form field to the group this button is part of.
 */
function addField(e) {
  // Find our parent field-group div
  let fieldGroup = this.parentElement

  // Get its keypath
  let keypath = fieldGroup.dataset.keypath

  // Find its last field child
  let existingFields = fieldGroup.getElementsByClassName('field')
  let templateField = existingFields[existingFields.length - 1]

  // Get its number
  let fieldNumber = templateField.dataset.number

  // Duplicate it
  let newField = templateField.cloneNode(true)

  // Increment the number and use the keypath and number to set IDs and cross
  // references.
  // TODO: Heavily dependent on the form field HTML. Maybe we want custom
  // elements for the labeled controls that know how to be list items?
  fieldNumber++
  newField.dataset.number = fieldNumber
  let newID = keypath + '[' + fieldNumber + ']'
  let newControl = newField.getElementsByClassName('control')[0]
  newControl.id = newID
  newControl.setAttribute('name', newID)
  let newLabel = newField.getElementsByTagName('label')[0]
  newLabel.setAttribute('for', newID)

  // Find the minus button
  let minusButton = fieldGroup.getElementsByClassName('remove-field')[0]

  // Put new field as a child before the minus button
  fieldGroup.insertBefore(newField, minusButton)

  // Enable the minus button
  minusButton.classList.remove('invisible')
}

/**
 * Remove the last form field from the group button is part of.
 */
function removeField(e) {
  // Find our parent field-group div
  let fieldGroup = this.parentElement

  // Find its field children
  let existingFields = fieldGroup.getElementsByClassName('field')

  if (existingFields.length > 1) {
    // There is a last field we can safely remove.
    let lastField = existingFields[existingFields.length - 1]
    fieldGroup.removeChild(lastField)
  }

  if (existingFields.length <= 1) {
    // Collection auto-updates. Now there's only one element. Don't let the
    // user remove it. If they don't want it, they can leave it empty.
    this.classList.add('invisible')
  }
}


// Change the submit button after hitting
function on_submit_button() {
    var f = document.getElementsByTagName('form')[0];
    if(f.checkValidity()) {
        var elem = document.getElementById("submit");
        elem.value = "Submitting...";
        elem.disabled = true;
        f.submit();
    } else {
        alert(document.getElementById('example').validationMessage);
        return false;
    }
}
