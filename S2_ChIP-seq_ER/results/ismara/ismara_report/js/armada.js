/*global $ */

function armada_get_remote_dir() {
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

function armada_send_data(remote_dir) {
    /* Put entered data into json format and sends it server */
    "use strict";
    console.log("send_data");
    var data = [];
    var rows =  $("#sample_index_table tbody tr");
    var rows_length = rows.length;
    var rpl_flag  = null;
    var replicates_map = {};
    var replicate_count = {};
    var i;
    for (i = 0; i < rows_length; i += 1) {
        var smp_name, smp_idx, rpl_idx, smp_new_name;
        smp_idx = $(rows[i]).find(".condition").val();
        /* all samples should be assigned to new datapoints */
        if (smp_idx === "0") {
            alert("Please set timepoints for all samples.");
            $(".armada_submit").show();
            return false;
        }
        console.log($(rows[i]).find("td a").html());
        smp_name = $(rows[i]).find("td a").html().replace(/^\s+|\s+$/g, "");
        console.log($(rows[i]).find(".timepoint_name input").val());
        smp_new_name = $(rows[i]).find(".timepoint_name input").val().replace(/\s+/g, '_');
        console.log($(rows[i]).find(".gap input").val());
        var gap = $(rows[i]).find(".gap input").val();
        if (gap == undefined) {
            gap = 0;
        }
        console.log("point1");
        console.log(gap);
        console.log($(rows[i]).find(".replicate").val());
        rpl_idx = $(rows[i]).find(".replicate").val();
        /* Check that replica indices are given either for all samples or for none */
        if (rpl_flag === null) {
            rpl_flag = Boolean(Number(rpl_idx));
        } else if (rpl_flag !== Boolean(Number(rpl_idx))) {
            alert("Please specify replicate indices either for all samples or for none.");
            $(".armada_submit").show();
            return false;
        }
        /* Check for repeating replicate number in one condiiton */
        if (replicates_map[smp_idx] === undefined) {
            replicates_map[smp_idx] = [rpl_idx];
        } else if (replicates_map[smp_idx].indexOf(rpl_idx) > -1 && rpl_idx > 0) {
            alert("You have given the same batch/replicate number more than once for one condiiton");
            $(".armada_submit").show();
            return false;
        } else {
            replicates_map[smp_idx].push(rpl_idx);
        }
        
        if (replicate_count[rpl_idx] === undefined) {
            replicate_count[rpl_idx] = 1;
        } else {
            replicate_count[rpl_idx] += 1;
        }

        /* Add data to list */
        data.push({"name": smp_name, "si": smp_idx, "ri": rpl_idx, "new_name": smp_new_name, "gap": gap});
    }

    for (var rpl_idx in replicate_count) {
        if (replicate_count.hasOwnProperty(rpl_idx)) {
            if (replicate_count[rpl_idx] === 1) {
                alert("It should be at least two batch/replicates with the same number in dataset. Perhaps you should try to submit your data without replicate indices.");
                $(".armada_submit").show();
                return false;
            }
        }
    }

    var post_url = "/fcgi/mara/armada";
    var secure = 0;
    if (document.URL.match(/\/secure\//)) {
        post_url = "/secure/fcgi/mara/armada";
        secure = 1;
    }
    $.post(post_url, {"dir": armada_get_remote_dir(),
                      "data": JSON.stringify(data),
                      "secure": secure
                     },
           function(x){
               $(".armada_submit").show();
               window.location = "armada";
           });
}


function armada_add_options(obj, type) {
    /* Add new option to select menu if last one was choosen */
    "use strict";
    var cur_val, num, child_num, i, class_selector;
    console.log($(obj).val(), obj, type);
    cur_val = parseInt($(obj).val(), 10);
    num = 0;
    child_num = $(obj).children().length;
    for (i = 0; i < child_num; i = i + 1) {
        var child = $($(obj).children()[i]);
        var value = parseInt(child.val(), 10);
        if (num < value) {
            num = value;
        }
    }
    class_selector = '.' + type;
    console.log("add_options: class_selector " + class_selector);
    if (num === cur_val) {
        num = num + 1;
        var label = type + num;
        if (type === 'replicate') {
            label = "replicate/batch" + num;
        }
        var text = "<option value=\"" + num + "\">" + label + "</option>";
        $(text).appendTo(class_selector);
    }
}

function armada_change_all(obj) {
    /* Cahnge name in all fields for the same sample */
    "use strict";
    var input = $(obj).children('input');
    var val = input.val();
    var mclass = input.attr("class");
    $("." + mclass).each(function(){$(this).val(val);});
}

function armada_change_gap(obj) {
    /* change gap values for corresponding timepoints */
    var input = $(obj).children("input");
    var val = input.val();
    var gclass = input.attr("class").match("gap_\\S+");
    console.log(gclass + " gclass")
    $("." + gclass).each(function(){$(this).val(val);});
}

function armada_change_timepoint_name(obj) {
    /* change condition_name values for corresponding timepoints */
    var input = $(obj).children("input");
    var val = input.val();
    var fclass = input.attr("class").match("name_\\S+");
    $("." + fclass).each(function(){$(this).val(val);});
}

function armada_add_name_input(obj, type) {
    /* adds text input for sample name */
    "use strict";
    var cur_val = $(obj).val();
    var name_field = $(obj).parent().children(".timepoint_name");
    var gap_field = $(obj).parent().children(".gap");
    console.log("add_name_input: name_field ");
    console.log(obj);
    var class_name = "name_" + type + cur_val;
    var cur_name = type + cur_val;
    if ($("." + class_name).length > 0) {
        cur_name = $("." + class_name).val();
    }
    var name_string = " Name " + "<input type=\"text\"  class=\"" + class_name +" fname\" name=\"fname\" value=\"" + cur_name +"\"/> ";
    name_field.html(name_string);
    if (class_name != "name_timepoint1") {
        console.log("add_name_input: class name " + class_name);
        var gap_class = "gap_" +  type + cur_val;
        var gap_string = "Time gap <input type=\"text\" class=\"" + gap_class +" gap\" name=\"gname\" value=\"1\" type=\"number\">";
        gap_field.html(gap_string);
    }
}

function armada_UrlExists(url)
{
    var http = new XMLHttpRequest();
    http.open('HEAD', url, false);
    http.send();
    return http.status!=404;
}

$(function () {
    
    var remote_dir = armada_get_remote_dir();
    if (remote_dir) {
        $("#armada_init").show();
        $("#armada_init").button();
    }
    /* Check if there is existing armada results */
    var url = "armada/index.html";
    $.ajax({
        url: url,
        success: function()
        {
            $("<div class=\"center\" style=\"font-size: 1.2em;text-decoration: underline;\"><a href=\"armada/index.html\">Armada results</a></div>").insertAfter("#armada_init");
        }
    });
    
    $("#armada_init").click(function() {
        $(".average").hide();
        $(".armada").show();
        $("th.average").html("Timepoint");
        /*$("select.condition").toggleClass("condition timepoint");
        $("select.timepoint").attr("name", "timrpoint")*/
        $(".armada_submit").show();
        $(".armada_submit").button();
        $(".armada_submit").click(function() {$(".armada_submit").hide();armada_send_data();});
        $(".condition").unbind();
        $(".condition").change(function() {armada_add_options(this, "timepoint");
                                           armada_add_name_input(this, "timepoint")});
        $("select.condition").html("<option value=\"0\">Select</option><option value=\"1\">timepoint1</option>");
        $(".replicate").change(function() {armada_add_options(this, "replicate");});
        $(".timepoint_name").keyup(function(){console.log("timepoint_nameg");armada_change_timepoint_name(this);})
        $(".gap").keyup(function(){console.log("gap");armada_change_gap(this)});
        window.location.hash = "sample_index";
    }); 
});
