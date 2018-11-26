/*global $ */

function search_gene_get_remote_dir() {
    /* returns remote results dir from url parsing */
    "use strict";
    var remote_dir, url, url_match;
    remote_dir = null;
    url = document.URL;
    url_match = url.match(/^(http|https):\/\/w*\.?(mara|ismara|bc2-web78.bc2)\.unibas\.ch.*\/(\S+)\/i?s?mara_report/);
    if (url_match) {
        remote_dir = url_match[3];
    }
    return remote_dir;
}

function search_gene_send_data() {
    /* Put entered data into json format and sends it server */
    "use strict";
    var term = $("#search_gene input[name='term']").val().trim();
    var dir = $("#search_gene input[name='dir']").val()

    if (term == "") {
        alert("Please enter gene name!");
        return;
    }
    
    var post_url = "/fcgi/search_gene";
    var secure = 0;
    if (document.URL.match(/\/secure\//)) {
        post_url = "/secure/fcgi/search_gene";
        secure = 1;
    }

    $.ajax({
        type: "POST",
        url: post_url, 
        data: {
            "dir": dir,
            "term": term,
            "secure": secure
        },
        beforeSend: function() {$("body").addClass("loading");},
        complete: function(){$("body").removeClass("loading");}
    })
        .done(
            function(x){
                console.log("x:" + x);
                if (x == "Fail") {
                    alert("Term " + term + " not found.");
                }
                else {
                    var scratch_url = "https://" + window.location.hostname;
                    /*scratch_url += "/ISMARA/scratch/"*/
                    var data = window.location.pathname.split("/");
                    scratch_url += data.slice(0,data.findIndex(function(x){return x == "ismara_report"}) - 1).join("/") + "/";
                    console.log(scratch_url + dir + "/" + term + ".html");
                    window.location = scratch_url + dir + "/" + term + ".html";
                }
            });
}

$(function () {
    var remote_dir = search_gene_get_remote_dir();
    if (remote_dir) {
        $("#search_gene").show();
//        $("#search_gene_button").button();
    }
    $("#search_gene_button").click(function(){search_gene_send_data();});
});
