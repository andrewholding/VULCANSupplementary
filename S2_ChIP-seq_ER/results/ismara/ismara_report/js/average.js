/*global $ */

function average_get_remote_dir() {
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

function average_send_data(remote_dir) {
    /* Put entered data into json format and sends it server */
    "use strict";
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
            alert("Please asign all samples to new condiitons.");
            return false;
        }
        smp_name = $(rows[i]).find("td a").html().replace(/^\s+|\s+$/g, "");
        smp_new_name = $(rows[i]).find(".condition_name input").val().replace(/\s+/g, '_');
        rpl_idx = $(rows[i]).find(".replicate").val();
        /* Check that replica indices are given either for all samples or for none */
        if (rpl_flag === null) {
            rpl_flag = Boolean(Number(rpl_idx));
        } else if (rpl_flag !== Boolean(Number(rpl_idx))) {
            alert("Please specify replicate indices either for all samples or for none.");
            return false;
        }
        /* Check for repeating replicate number in one condiiton */
        if (replicates_map[smp_idx] === undefined) {
            replicates_map[smp_idx] = [rpl_idx];
        } else if (replicates_map[smp_idx].indexOf(rpl_idx) > -1 && rpl_idx > 0) {
            alert("You have given the same batch/replicate number more than once for one condiiton");
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
        data.push({"name": smp_name, "si": smp_idx, "ri": rpl_idx, "new_name": smp_new_name});
    }

    for (var rpl_idx in replicate_count) {
        if (replicate_count.hasOwnProperty(rpl_idx)) {
            if (replicate_count[rpl_idx] === 1) {
                alert("It should be at least two batch/replicates with the same number in dataset. Perhaps you should try to submit your data without replicate indices.");
                return false;
            }
        }
    }

    var post_url = "/fcgi/mara/average";
    var secure = 0;
    if (document.URL.match(/\/secure\//)) {
        post_url = "/secure/fcgi/mara/average";
        secure = 1;
    }
    $.post(post_url, {
        "dir": average_get_remote_dir(),
        "data": JSON.stringify(data),
        "email": $("#email").val(),
        "project_name": $("#project_name").val(),
        "orgnl_link": $("#orgnl_link").val(),
        "secure": secure
    },
           function(x){
               var scratch_url = "https://" + window.location.hostname 
               if (secure === 1) {
                   scratch_url += "/secure"
               }
               scratch_url += "/ISMARA/scratch/"
               window.location = scratch_url + x + "/averaged_report";
                      });
}


function average_add_options(obj, type) {
    /* Add new option to select menu if last one was choosen */
    "use strict";
    var cur_val, num, child_num, i, class_selector;
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

function average_change_all(obj) {
    /* Cahnge name in all fields for the same sample */
    "use strict";
    var input = $(obj).children('input');
    var val = input.val();
    var mclass = input.attr("class");
    $("." + mclass).each(function(){$(this).val(val);});
}

function average_add_name_input(obj, type) {
    /* adds text input for sample name */
    "use strict";
    var cur_val = $(obj).val();
    var name_field = $(obj).parent().children(".condition_name");
    var class_name = "name_" + type + cur_val;
    var cur_name = type + cur_val;
    if ($("." + class_name).length > 0) {
        cur_name = $("." + class_name).val();
    }
    name_field.html(" Name " + "<input type=\"text\"  class=\"" + class_name +"\" name=\"fname\" value=\"" + cur_name +"\"/>");
}

$(function () {
    /* averaging */
    var remote_dir = average_get_remote_dir();
    if (remote_dir) {
        $("#average_init").show();
        $("#average_init").button();
    }

    $("#average_init").click(function(){
        $(".armada").hide();
        $(".average").show();
        $(".chipseq-tools").hide();
        $(".average_submit").button();
        $(".average_submit").click(function() { average_send_data();});
        $("#batch_button").button();
        $("#batch_button").click(function(){$(".batch").show();});
        $(".condition").unbind();
        $(".condition").change(function() {average_add_options(this, "condition");
                                           average_add_name_input(this, "condition")});
        $(".replicate").change(function() {average_add_options(this, "replicate");});
        $(".condition_name").keyup(function(){average_change_all(this);})
        window.location.hash = "sample_index";
    });
    $(".average_show_conf").button();
    $(".average_hide_conf").button();
    $("#avrg_conf").hide();
    $(".average_show_conf").click(function(){$("#avrg_conf").show();
                                             $(".average_hide_conf").show();
                                             $(".average_show_conf").hide();
                                            });
    $(".average_hide_conf").click(function(){$("#avrg_conf").hide();
                                             $(".average_hide_conf").hide();
                                             $(".average_show_conf").show();
                                             window.location.hash = "average_configuration";
                                            });
       
    
});
