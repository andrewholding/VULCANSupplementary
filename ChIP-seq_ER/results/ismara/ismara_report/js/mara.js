/*global $:false */

function parseURL(url) {

	var parser = document.createElement('a'),
		searchObject = {},
		queries, split, i;

	// Let the browser do the work
	parser.href = url;

	// Convert query string to object
	queries = parser.search.replace(/^\?/, '').split('&');
	for( i = 0; i < queries.length; i++ ) {
		split = queries[i].split('=');
		searchObject[split[0]] = split[1];
	}

	return {
		protocol: parser.protocol,
		host: parser.host,
		hostname: parser.hostname,
		port: parser.port,
		pathname: parser.pathname,
		search: parser.search,
		searchObject: searchObject,
		hash: parser.hash
	};
}

function get_protocol(data_url) {
    var parser = parseURL(data_url);
    return parser.protocol;
}

function get_filename(data_url) {
    var parser = parseURL(data_url);
    return parser.pathname.split("/").pop();
}

function validate_user_data() {
    var validate_result = true;

    if ($("#type_rnaseq").attr("checked") === "checked" || $("#type_chipseq").attr("checked") === "checked") {
        if ($("#organism_row input:checked").length === 0) {
            alert("Please choose one of organisms.");
            $("#organism_row").css("border", "2px solid red");
            validate_result = false;
        }
    }
    var data_url = $("#url").val();
    if (data_url === undefined) {
        data_url = "";
    }
    if (data_url.trim() !== "") {
        var filename = get_filename(data_url);
        if (!filename.match(/(tar\.gz|\.tar|\.zip|\.tgz)$/i)) {
            alert("Could not find valid filename from your URL. It should be .tar, .tgz, .tar.gz or .zip file.");
            validate_result = false;
        }
    }
    return validate_result;
}

function submit_user_data() {
    var email = $("#email").val();
    var project = $("#project").val();
    var type = $("#type_row input:checked").val();
    var organism = $("#organism_row input:checked").val() === undefined ? "" : $("#organism_row input:checked").val();
    var mirna = $("#mirna_row input:checked").val();
    var data_url = $("#url").val();
    var protocol = "";
    if (data_url === undefined || data_url.trim() === "") {
        data_url = "";
    } else {
        protocol = data_url.split("://")[0];
        data_url = data_url.split("://")[1];
    }

    var submit_url = location.pathname + "/run"
    if ($("#logged").val() === "1") {
        submit_url = "/secure" + submit_url
    }

    $(".container").css('height', $(".container").css('height'));
    $(".container").css('width', $(".container").css('width'));
    $(".container").css('border', 'none');
    $(".container").html("<img src=\"/img/loader.gif\" alt=\"loader\">");

    img_margin_top = Math.round($(".container").height()/2 - 33);
    img_margin_left =  Math.round($(".container").width()/2 - 33);
    $(".container img").css('margin-top', img_margin_top + "px");
    $(".container img").css('margin-left', img_margin_left + "px");

    $.post(submit_url, {
        "email": email,
        "project": project,
        "type": type,
        "organism": organism,
        "mirna": mirna,
        "data_url": data_url,
        "protcol": protocol
    },
           function(x){
               if(x.match(/data_\S+$/) !== null) {
                   window.location = "https://" + window.location.hostname + x;
               }
               else {
                   $(".container").html(x+"<br><br><h2><a href=''>Go back to main page</a></h2>");
               }
           });
    
}

function validateExpertForm() {
    /* Test if at least one file is submitted */
    if ($("#expr_table").val() === "" || $("#sitecount_table").val() === "") {
        alert("Please provide one file with expression table and one file with sitecount matrix.");
        return false;
    }
    return true;
}
    
function change_fileupload_url() {
    var curr_url = $("#fileupload").attr("action");
    var new_url = curr_url.replace("/fcgi/", "/secure/fcgi/");
    $("#fileupload").attr("action", new_url);
}

function login() {
    var username = $("#username").val();
    var password = $("#password").val();
    var posting = $.post("/secure" + location.pathname + "/login", {
        "username": username,
        "password": password
    }, function(x) {
        if (x === "Log in failed.") {
            $("#login #status").html(x);
        }
        else {
            $("#in").hide();
            $("#login #logout").show();
            var data = jQuery.parseJSON(x);
            $("#login #status").html("Welcome " + data[0].username + "<ul><li><a href=\"/secure" + location.pathname + "/list\">List jobs</a></li><li><a href=\"" + location.pathname + "/pchng\">Change password</a></li></ul>");
            $("#email").val(data[0].email);
            $("#logged").val("1");
            change_fileupload_url();
        }        

    });
}

function logout() {
    var posting = $.post(location.pathname + "/logout",
                         {},
                         function(x) {
                             if (x === "Ok") {
                                 $("#in").show();
                                 $("#login #logout").hide();
                                 $("#login #status").html("");
                                 $("#logged").val("0");
                             } else {
                                 $("#login #status").html("Failed to log out. Please clean the coockies and reload the page.");
                             }
                         });
}

$(function () {
    'use strict';
    if (window.location.protocol === "http:") {
        window.location = "https://" + window.location.hostname + window.location.pathname;

    }
    var flag = $("#flag").val();
    if (flag === "1") window.location.reload();
    $("#flag").val("1");
    var files4upload = [];
    var uploaded_files = 0;

    /* Toggling TOC items */
    $(".toggle").click(function(){
        var info_block_id  = $(this).attr("id").split("_");
        $("#info div").hide();
        $(".toggle").css("text-decoration", "none");
        $("#" + info_block_id).show();
        $(this).css("text-decoration", "underline");
        
    });

    $("#type_microarray").click(function(){
        $("#organism_row").hide();
        $("#mirna_true").attr("checked", true);
        $("#mirna_false").attr("checked", false);
        $("#mirna_row").show();
        $("#fileupload").show();
        $("#expert_form").hide();
        $("#fastq_help").hide();
    });
    $("#type_rnaseq").click(function(){
        $("#organism_row").show();
        $("#mirna_true").attr("checked", true);
        $("#mirna_false").attr("checked", false);    
        $("#mirna_row").show();
        $("#fileupload").show();
        $("#expert_form").hide();
        $("#fastq_help").show();
        if (window.location.href.search(/\?mm10lncrna/) != -1) {
            $(".mm10lncrna").show();
        }
    });
    $(".mm10lncrna input").click(function(){
            $("#mirna_true").attr("checked", false);
            $("#mirna_false").attr("checked", true);    
            $("#mirna_row").hide();                    
    });

    $("#type_chipseq").click(function(){
        $("#organism_row").show();
        $("#mirna_true").attr("checked", false);
        $("#mirna_false").attr("checked", true);
        $("#mirna_row").show();
        $("#fileupload").show();
        $("#expert_form").hide();
        $("#fastq_help").show();
    });
    $("#type_expert").click(function(){
        $("#organism_row").hide();
        $("#mirna_row").hide();
        $("#fileupload").hide();
        $("#expert_form").show();
        $("#fastq_help").hide();
    });
    $("#submit_user_data").click(function(){
        var validate_result = validate_user_data();
        if (validate_result) {
            submit_user_data();
        }
    });

    $("#urlSubmit").click(function(){
        var validate_result = validate_user_data();
        if (validate_result) {
            submit_user_data();
        }
    });
    $("#login #authenticate").click(function(){
        login();
    });
    $("#login #logout").click(function(){
        logout();
    });

    /* fileupload bindings */
    $('#fileupload').bind('fileuploadadd', function (e, data) {
        if (files4upload.indexOf(data.files[0].name) !== -1){
            data.files[0].error = "Duplicated";
        }
        var post_url = location.pathname + "/swd";
        if ($("#logged").val() === "1") {
            post_url = "/secure" + post_url;
        }
        var swd = $("#swd").val();
        if (swd === "0") {
            $.get(post_url, function(data){
                $("#swd").val("1");
            });
        }
    });
    $('#fileupload').bind('fileuploadadded', function (e, data) {
        if (data.files[0].error === null) {
            files4upload.push(data.files[0].name);
        }
    });

    $('#fileupload').bind('fileuploadfailed', function (e, data) {
        if (data.files[0].error === null) {
            files4upload.splice(files4upload.indexOf(data.files[0].name), 1);
        }
    });
    
    $('#fileupload').bind('fileuploadstarted', function (e, data) {
        $(".cancel").hide();
    });

    $('#fileupload').bind('fileuploadcompleted', function (e, data) {
        uploaded_files += 1;
        if (uploaded_files === files4upload.length) {
            var validate_result = validate_user_data();
            if (validate_result) {
                submit_user_data();
            }
            else {
                $(".span7").hide();
                $("#submit_user_data").show();
                $("#brief_instructions").html("Please fill missing values and click \"Submit\" button.");
            }
        }
    });

    $('#fileupload').fileupload({
        maxChunkSize: 50000000,
        sequentialUploads: true,
        acceptFileTypes: /(\.|\/)(bed|cel|gz|bz2|zip|tar|bam|sam|fastq)$/i
    });
    $(".tip").tooltip({showURL: false, track: true});
});
