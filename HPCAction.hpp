//------------------------------------------------------------------------------
/// @file
/// @brief    Action 列挙型
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  Copyright (c) 2015 HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください

//------------------------------------------------------------------------------
#pragma once

namespace hpc {

    //------------------------------------------------------------------------------
    /// @brief 行動を表します。
    enum Action {
        Action_MoveLeft,  ///< 左に移動。x - 1。
        Action_MoveRight, ///< 右に移動。x + 1。
        Action_MoveDown,  ///< 下に移動。y - 1。
        Action_MoveUp,    ///< 上に移動。y + 1。

        Action_TERM
    };
}
//------------------------------------------------------------------------------
// EOF