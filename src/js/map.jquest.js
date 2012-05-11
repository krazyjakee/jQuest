var Map = {
	settings: {
		fileExtension: 'json'
	},
	mapData: false,
	drawMap: function(targetElem){
		$.each(Map.mapData.layers, function(index, layer){
			if(layer.type == "tilelayer"){
				$(targetElem).append('<div class="layer" id="layer'+index+'"></div>');
				for(var i = 0; i < layer.data.length; i++){
					$('#layer'+index).append('<div id="tile'+index+'-'+i+'" class="tile"></div>');
					Map.drawTile(index, layer.data[i], i);
				}
				$('#layer'+index).css('width',Map.mapData.tilewidth * Map.mapData.width + 'px');
				$('#layer'+index).css('height',Map.mapData.tileheight * Map.mapData.height + 'px');
				$(targetElem).addClass('viewport');
				$('#maploading').remove();
			}else if(layer.type == "objectlayer"){
				
			}
		});
	},
	loadMap: function(targetElem, mapSource, callback){
		if(typeof mapSource == 'object'){
			Map.mapData = mapSource;
			Map.drawMap(targetElem);
			if(callback){
				callback();
			}
		}else if(typeof mapSource == 'string'){
			$.getJSON(mapSource+'.'+Map.settings.fileExtension, function(json){
				Map.mapData = json;
				Map.drawMap(targetElem);
				if(callback){
					callback();
				}
			});
		}
	},
	drawTile: function(layer, srcTile, targetTile){
		var target = $('#tile'+layer+'-'+targetTile);
			
		var setData = false;
		$.each(Map.mapData.tilesets, function(index, tileset){
			if(srcTile >= tileset.firstgid){
				setData = tileset;
			}else{
				return false;
			}
		});
		
		var setWidth = Math.floor(setData.imagewidth / setData.tilewidth)-1;
		var setHeight = Math.floor(setData.imageheight / setData.tileheight)-1;
		
		var heightCount = 0;
		var widthCount = 0;
		for(var i = 1; i < srcTile; i++){
			if(widthCount < setWidth){
				widthCount++;
			}else{
				widthCount = 0;
				heightCount++;
			}
		}
		
		var offset = {
			x: (widthCount * setData.tilewidth) + (setData.spacing * widthCount) + setData.margin,
			y: (heightCount * setData.tileheight) + (setData.spacing * heightCount) + setData.margin
		}
		
		$(target).width(setData.tilewidth);
		$(target).height(setData.tileheight);
		$(target).css('background-image','url('+setData.image+')');
		$(target).css('background-position','-'+offset.x+'px -'+offset.y+'px');
	}
}